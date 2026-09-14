% A parallel wrapper around GRAPE that enables ensemble optimal control
% optimisations. This function handles systems with multiple control po-
% wer levels, multiple resonance offsets, multistate transfers, ensemb-
% les of drift Liouvillians, etc. Syntax:
%
%          [traj_data,fidelity,...
%           gradient,hessian]=ensemble(waveform,spin_system)
%
% Parameters:
%
%   waveform  - control coefficients for each control operator, rad/s
%
% Outputs:
%
%   traj_data    - trajectory data for subsequent diagnostic plotting
%
%   fidelity     - figure of merit for the overlap of the current state
%                  of the system and the desired state(s). When penalty
%                  methods are specified, fidelity is returned as an ar-
%                  ray separating the penalties from the simulation
%                  fidelity.
%
%   gradient     - gradient of the fidelity with respect to the control
%                  sequence. When penalty methods are specified, gradi-
%                  ent is returned as an array separating penalty gra-
%                  dients from the fidelity gradient.
%
%   hessian      - Hessian of the fidelity with respect to the control
%                  sequence. When penalty methods are specified, gradi-
%                  ent is returned as an array separating penalty Hes-
%                  sians from the fidelity Hessian.
%
% Note: the ensemble cases enumerated in spin_system.control.catalog
%       are processed in the contiguous per-worker blocks assigned by
%       optimcon.m in spin_system.control.worker_cases. Each worker
%       holds the common frozen problem and the drift generators of its
%       own block, published by optimcon.m as pool constants, and grafts
%       the live client-side control structure on top of them, so only
%       the waveform and the live control fields travel at each objec-
%       tive evaluation; the gradient, the Hessian, and averaged trajec-
%       tories are summed on the workers. This func-
%       tion must be called from the client, on the pool that was
%       open when optimcon.m ran: a worker holds only its own block.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ensemble.m>

function [traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)

% Check consistency
grumble(spin_system,waveform);

% Worker-resident problem data handles
invariants=spin_system.control.invariants;
drift_slices=spin_system.control.drift_slices;

% Live problem data is the client-side control structure less what the workers already hold
control=rmfield(spin_system.control,intersect({'invariants','drift_slices','worker_cases','basis'},...
                                               fieldnames(spin_system.control)));
control.return_traj=isfield(control,'return_traj')&&control.return_traj;

% Count the outputs and the cases
n_outputs=nargout; n_cases=size(control.catalog,1);
if (n_outputs>3)&&(~all(cellfun(@(f)isequal(f,@no_dist),control.distortion(:))))
    error('Hessians are not available with waveform distortions.');
end

% Run the ensemble loop, each worker over its own case block
spmd (poolsize)

    % Evaluate the block of cases assigned to this worker
    [traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value,...
                                                           control,spmdIndex,waveform,n_outputs);

    % Reduce to the first worker and pack
    results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1),...
                   'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1));

end

% Collect from the first worker
results=results{1}; traj_data=results.traj; fidelities=results.fid;
gradient=results.grad; hessian=results.hess;

% Average the block trajectory sums
if ismember('average',control.traj_opts)
    ave_traj=traj_data{1}.forward;
    for n=2:numel(traj_data)
        ave_traj=ave_traj+traj_data{n}.forward;
    end
    traj_data={struct('forward',{(1/n_cases)*ave_traj})};
end

% Ensemble averages of fidelity, gradient, and Hessian
fidelity=sum(fidelities)/n_cases;
if n_outputs>2
    gradient=reshape(gradient/n_cases,size(waveform));
end
if (n_outputs>3)&&strcmp(control.integrator,'rectangle')
    hessian=reshape(hessian/n_cases,numel(waveform)*[1 1]);
else
    hessian=[];
end

% Run diagnostic plotting (expensive!)
if ~isempty(spin_system.control.plotting)

    % With or without instrumental distortions
    if ~isempty(spin_system.control.distplot)

        % Apply the distortions
        dist_waveform=waveform;
        for k=1:numel(spin_system.control.distplot)

            % Extract and apply distortion function
            dist_function=spin_system.control.distplot{k};
            dist_waveform=dist_function(dist_waveform);

        end

        % Real-life trajectory and the distorted control sequence
        ctrl_trajan(spin_system,dist_waveform,traj_data,fidelities);

    else

        % Real-life trajectory but the ideal control sequence
        ctrl_trajan(spin_system,waveform,traj_data,fidelities);

    end

end

end

% Fidelity, gradient, and Hessian contributions of one block of ensemble cases
function [traj,fid,grad,hess]=ens_block(ss,drifts,control,block,waveform,n_outputs)

% Graft live client data over the frozen worker copy, keep what only the worker holds
frozen=ss.control; ss.control=control; ss.control.drifts=drifts;
missing=setdiff(fieldnames(frozen),fieldnames(control));
for k=1:numel(missing)
    ss.control.(missing{k})=frozen.(missing{k});
end

% GRAPE function for the formalism
switch ss.bas.formalism
    case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}
        grape=@grape_liouv;
    case 'zeeman-hilb'
        grape=@grape_hilb;
    otherwise
        error('unrecognised formalism specification.');
end

% Cases of this block, waveform dimensions, and offset ensemble size
my_cases=frozen.worker_cases{block}; n_mine=numel(my_cases); catalog=control.catalog;
ncont=size(waveform,1); nsteps=size(waveform,2); off_ens_sizes=cellfun(@numel,control.offsets);

% Preallocate block outputs, derivative buffers only when requested
traj=cell(n_mine,1); fid=zeros(1,n_mine); grad=[]; hess=[];
if n_outputs>2, grad=zeros(ncont*nsteps,1); end
if n_outputs>3, hess=zeros((ncont*nsteps)^2,1); end

% Loop over the cases of the block
for m=1:n_mine

    % Extract ensemble indices
    n=my_cases(m); n_rho=catalog(n,1); n_sys=catalog(n,2);
    n_pwr=catalog(n,3); n_off=catalog(n,4);
    n_phi=catalog(n,5); n_dis=catalog(n,6);

    % Get initial and target states, drift, and waveform
    rho_init=control.rho_init{n_rho}; rho_targ=control.rho_targ{n_rho};
    L=ss.control.drifts{n_sys}; local_waveform=waveform;

    % Phase cycle: a rotation of each control channel and phases on the states
    R=eye(ncont);
    if ~isempty(control.phase_cycle)
        phi=control.phase_cycle(n_phi,:);
        rho_init=exp(1i*phi(1))*rho_init; rho_targ=exp(1i*phi(end))*rho_targ;
        R=kron(diag(cos(phi(2:end-1))),eye(2))+kron(diag(sin(phi(2:end-1))),[0 -1; 1 0]);
        local_waveform=R*local_waveform;
    end

    % Add offset terms, first channel index fastest (user specifies offsets in Hz)
    if ~isempty(off_ens_sizes)
        off_idx=cell(1,numel(off_ens_sizes)); [off_idx{:}]=ind2sub([off_ens_sizes 1],n_off);
        for k=1:numel(off_ens_sizes)
            L=L+sparse(2*pi*control.offsets{k}(off_idx{k})*ss.control.off_ops{k});
        end
    end

    % Move the waveform into physical units
    power_lvl=control.pwr_levels(n_pwr); local_waveform=power_lvl*local_waveform;

    % Apply waveform distortions, with their Jacobian when derivatives are needed
    if n_outputs>2, J=speye(numel(local_waveform)); end
    for k=1:size(control.distortion,2)
        if n_outputs>2
            [local_waveform,stage_jacobian]=control.distortion{n_dis,k}(local_waveform);
            J=stage_jacobian*J;
        else
            local_waveform=control.distortion{n_dis,k}(local_waveform);
        end
    end

    % Fidelity, trajectory, and derivatives
    outputs=cell(1,n_outputs);
    [outputs{:}]=grape(ss,L,ss.control.operators,local_waveform,rho_init,rho_targ,control.fidelity);
    traj{m}=outputs{1}; fid(m)=outputs{2};

    % Gradient through the Jacobian, the phase cycle, and the power level
    if n_outputs>2
        grad_n=R'*reshape(J'*outputs{3}(:),ncont,nsteps);
        grad=grad+power_lvl*grad_n(:);
    end

    % Hessian through the phase cycle and the power level
    if n_outputs>3
        K=kron(speye(nsteps),sparse(R')); hess_n=K*outputs{4}*K';
        hess=hess+power_lvl^2*hess_n(:);
    end

end

% Collapse a non-empty block into one trajectory sum when only the average is needed
if ismember('average',control.traj_opts)&&(n_mine>0)
    traj_sum=traj{1}.forward;
    for m=2:n_mine
        traj_sum=traj_sum+traj{m}.forward;
    end
    traj={struct('forward',{traj_sum})};
end

end

% Consistency enforcement
function grumble(spin_system,waveform)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if ~all(isfield(spin_system.control,{'catalog','ens_sizes','invariants','drift_slices','frozen_fields','worker_cases','pool_id'}))
    error('ensemble catalog missing from spin_system, run optimcon() first.');
end
if ~isempty(getCurrentWorker())
    error('ensemble() must run on the client: the frozen problem is distributed over the pool workers.');
end
current_pool=gcp('nocreate'); pool_id=0;
if ~isempty(current_pool), pool_id=current_pool.ID; end
if pool_id~=spin_system.control.pool_id
    error('parallel pool changed after optimcon(), re-run optimcon().');
end
if any(isfield(spin_system.control,spin_system.control.frozen_fields))
    error('generators and operators are frozen after optimcon(), re-run optimcon() to change them.');
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be an array of real numbers.');
end
if size(waveform,1)~=spin_system.control.ncontrols
    error('the number of rows in waveform must equal to the number of controls.');
end
switch spin_system.control.integrator
    case 'rectangle'
        if size(waveform,2)~=spin_system.control.pulse_nsteps
            error('the number of columns in waveform must be equal to the number of time steps.');
        end
    case 'trapezium'
        if size(waveform,2)~=(spin_system.control.pulse_nsteps+1)
            error('the number of columns in waveform must be (number of time steps)+1.');
        end
    otherwise
        error('unknown time propagation algorithm.');
end
if numel(spin_system.control.pulse_dt)~=spin_system.control.pulse_nsteps
    error('the length of control.pulse_dt has changed, re-run optimcon().');
end
[catalog_now,ens_sizes_now]=ens_catalog(spin_system.control);
if (numel(spin_system.control.rho_targ)~=numel(spin_system.control.rho_init))||...
   (~isequal(ens_sizes_now,spin_system.control.ens_sizes))||...
   (~isequal(catalog_now,spin_system.control.catalog))
    error('ensemble composition changed after optimcon(), re-run optimcon().');
end
end

% "After unsuccessful attempts to trap redtail monkeys
% at the Zika Forest with the intention of live-bleeding
% and release, monkeys had to be sampled by means of
% 12-bore shotguns."
%
% https://doi.org/10.1016/0035-9203(82)90161-4

