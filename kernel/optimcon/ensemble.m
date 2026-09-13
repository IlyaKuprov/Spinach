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
%       holds the frozen problem data of its own block, published by
%       optimcon.m as a pool constant, and grafts the live client-side
%       control structure on top of it, so only the waveform and the
%       live control fields travel at each objective evaluation; the
%       gradient and the Hessian are summed on the workers. This func-
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

% Pull the worker-resident problem data handle
invariants=spin_system.control.invariants;

% Live problem data is the client-side control structure
control=rmfield(spin_system.control,'invariants');

% Default the trajectory return flag
control.return_traj=isfield(control,'return_traj')&&control.return_traj;

% Count the outputs and the cases
n_outputs=nargout; n_cases=size(control.catalog,1);

% Run the ensemble loop, each worker over its own case block
spmd (poolsize)

    % Evaluate the block of cases assigned to this worker
    [traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,control,...
                                                           control.worker_cases{spmdIndex},...
                                                           waveform,n_outputs);

    % Reduce to the first worker and pack
    results=struct('traj',{spmdCat(traj_local,1,1)},'fid',{spmdCat(fid_local,2,1)},...
                   'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1));

end

% Collect from the first worker
results=results{1}; traj_data=results.traj; fidelities=results.fid;
gradient=results.grad; hessian=results.hess;

% Apply trajectory options
if ismember('average',spin_system.control.traj_opts)

    % Average the block trajectory sums
    ave_traj=(1/n_cases)*traj_data{1}.forward;
    for n=2:numel(traj_data)
        ave_traj=ave_traj+(1/n_cases)*traj_data{n}.forward;
    end

    % Overwrite traj_data
    traj_data=[]; traj_data{1}.forward=ave_traj;

end

% Add up fidelities
fidelities=cell2mat(fidelities);
fidelity=sum(fidelities)/n_cases;

% Normalise gradient
if n_outputs>2
    gradient=reshape(gradient/n_cases,size(waveform));
end

% Normalise Hessian
if (n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')
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
function [traj,fid,grad,hess]=ens_block(ss,control,my_cases,waveform,n_outputs)

% Graft live client data over the frozen worker copy
frozen=ss.control; ss.control=control;
for k=1:numel(control.frozen_fields)
    fname=control.frozen_fields{k};
    ss.control.(fname)=frozen.(fname);
end

% Waveform dimension statistics
ncont=size(waveform,1); nsteps=size(waveform,2);

% Case catalog and offset ensemble size
catalog=control.catalog; off_ens_sizes=cellfun(@numel,control.offsets);

% Number of cases in the block
n_mine=numel(my_cases);

% Preallocate local outputs
traj=cell(n_mine,1); fid=cell(1,n_mine);
grad=zeros(ncont*nsteps,1); hess=[];
if n_outputs>3, hess=zeros((ncont*nsteps)^2,1); end

% Loop over the local cases
for m=1:n_mine

    % Extract ensemble indices
    n=my_cases(m); n_rho=catalog(n,1); n_sys=catalog(n,2);
    n_pwr=catalog(n,3); n_off=catalog(n,4);
    n_phi=catalog(n,5); n_dis=catalog(n,6);

    % Get initial and target state
    rho_init=control.rho_init{n_rho};
    rho_targ=control.rho_targ{n_rho};

    % Localise the waveform
    local_waveform=waveform;

    % Apply the phase cycle
    if ~isempty(control.phase_cycle)

        % Apply phase to the initial state
        phi=control.phase_cycle(n_phi,1);
        rho_init=exp(1i*phi)*rho_init;

        % Apply phase to the target state
        phi=control.phase_cycle(n_phi,end);
        rho_targ=exp(1i*phi)*rho_targ;

        % Apply phases to the waveform
        for k=1:(size(local_waveform,1)/2)

            % Assemble complex waveform
            cplx_wave=local_waveform(2*k-1,:)+...
                   1i*local_waveform(2*k,:);

            % Get the phase
            phi=control.phase_cycle(n_phi,k+1);

            % Apply the phase
            cplx_wave=exp(1i*phi)*cplx_wave;

            % Get back X and Y components
            local_waveform(2*k-1,:)=real(cplx_wave);
            local_waveform(2*k,:)=imag(cplx_wave);

        end

    end

    % Get the drift generators
    L=ss.control.drifts{n_sys};

    % Add offset terms
    if ~isempty(off_ens_sizes)

        % Multi-index mathematics
        cum_sizes=fliplr(cumprod(off_ens_sizes));
        cum_sizes=[cum_sizes(2:end) 1]; lin_idx=n_off;
        for k=1:numel(off_ens_sizes)

            % Current channel index
            vi=rem(lin_idx-1,cum_sizes(k))+1;
            vj=(lin_idx-vi)/cum_sizes(k)+1;
            oper_idx=numel(off_ens_sizes)-k+1;

            % Add to the drift (user specifies offsets in Hz)
            L=L+sparse(2*pi*control.offsets{oper_idx}(vj)*...
                            ss.control.off_ops{oper_idx});

            % Next channel
            lin_idx=vi;

        end

    end

    % Move the waveform into physical units
    power_lvl=control.pwr_levels(n_pwr);
    local_waveform=power_lvl*local_waveform;

    % Call GRAPE
    if n_outputs==2

        % Apply waveform distortions
        for k=1:size(control.distortion,2)

            % Get distortion function
            dist_function=control.distortion{n_dis,k};

            % Apply distortion function
            local_waveform=dist_function(local_waveform);

        end

        % Fidelity and trajectory
        switch ss.bas.formalism

            case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}

                % Call Liouville space version of the GRAPE function
                [traj{m},fid{m}]=grape_liouv(ss,L,ss.control.operators,...
                                                         local_waveform,rho_init,rho_targ,...
                                                         control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj{m},fid{m}]=grape_hilb(ss,L,ss.control.operators,...
                                                        local_waveform,rho_init,rho_targ,...
                                                        control.fidelity);

            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end

    elseif n_outputs==3

        % Get the Jacobian going
        J=speye(numel(local_waveform));

        % Apply waveform distortions
        for k=1:size(control.distortion,2)

            % Get distortion function
            dist_function=control.distortion{n_dis,k};

            % Apply distortion and get its Jacobian
            [local_waveform,stage_jacobian]=dist_function(local_waveform);

            % Combine Jacobians
            J=stage_jacobian*J;

        end

        % Fidelity and trajectory
        switch ss.bas.formalism

            case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}

                % Call Liouville space version of the GRAPE function
                [traj{m},fid{m},grad_n]=grape_liouv(ss,L,ss.control.operators,...
                                                                      local_waveform,rho_init,rho_targ,...
                                                                      control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj{m},fid{m},grad_n]=grape_hilb(ss,L,ss.control.operators,...
                                                                     local_waveform,rho_init,rho_targ,...
                                                                     control.fidelity);

            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end

        % Store the gradient layout
        [n_rows,n_cols]=size(grad_n);

        % Stretch and apply the Jacobian
        grad_n=J'*grad_n(:);

        % Restore the original gradient layout
        grad_n=reshape(grad_n,[n_rows n_cols]);

    elseif n_outputs==4

        % Fidelity and trajectory
        switch ss.bas.formalism

            case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}

                % Call Liouville space version of the GRAPE function
                [traj{m},fid{m},...
                 grad_n,hess_n]=grape_liouv(ss,L,ss.control.operators,...
                                                       local_waveform,rho_init,rho_targ,...
                                                       control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj{m},fid{m},...
                 grad_n,hess_n]=grape_hilb(ss,L,ss.control.operators,...
                                                      local_waveform,rho_init,rho_targ,...
                                                      control.fidelity);

            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end

    end

    % Post-process gradient
    if (~isempty(control.phase_cycle))&&(n_outputs>2)

        % Un-apply phases to gradient
        for k=1:(size(grad_n,1)/2)

            % Assemble complex gradient
            cplx_grad=grad_n(2*k-1,:)+...
                   1i*grad_n(2*k,:);

            % Get the phase
            phi=control.phase_cycle(n_phi,k+1);

            % Un-apply the phase
            cplx_grad=exp(-1i*phi)*cplx_grad;

            % Get back X and Y components
            grad_n(2*k-1,:)=real(cplx_grad);
            grad_n(2*k,:)=imag(cplx_grad);

        end

    end

    % Post-process Hessian
    if (~isempty(control.phase_cycle))&&(n_outputs>3)

        % Re-shape the Hessian as [ncont x nsteps x nsteps x ncont]
        hess_n=reshape(hess_n,[ncont nsteps ncont nsteps]);

        % Un-apply phases to Hessian
        for k=1:(size(grad_n,1)/2)

            % Get the phase
            phi=control.phase_cycle(n_phi,k+1);

            % Assemble complex Hessian - left
            cplx_hess=hess_n(2*k-1,:,:,:)+...
                   1i*hess_n(2*k,:,:,:);

            % Un-apply the phase
            cplx_hess=exp(-1i*phi)*cplx_hess;

            % Get back X and Y components
            hess_n(2*k-1,:,:,:)=real(cplx_hess);
            hess_n(2*k,:,:,:)=imag(cplx_hess);

            % Assemble complex Hessian - right
            cplx_hess=hess_n(:,:,2*k-1,:)+...
                   1i*hess_n(:,:,2*k,:);

            % Un-apply the phase
            cplx_hess=exp(-1i*phi)*cplx_hess;

            % Get back X and Y components
            hess_n(:,:,2*k-1,:)=real(cplx_hess);
            hess_n(:,:,2*k,:)=imag(cplx_hess);

        end

        % Reshape the Hessian back
        hess_n=reshape(hess_n,[ncont*nsteps nsteps*ncont]);

    end

    % Apply power level and accumulate
    if n_outputs>2
        grad=grad+power_lvl*grad_n(:);
    end
    if n_outputs>3
        hess=hess+power_lvl*power_lvl*hess_n(:);
    end

end

% Collapse the block into one trajectory sum when only the average is needed
if ismember('average',control.traj_opts)
    traj_sum=0;
    for m=1:n_mine
        traj_sum=traj_sum+traj{m}.forward;
    end
    traj={struct('forward',traj_sum)};
end

end

% Consistency enforcement
function grumble(spin_system,waveform)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if ~all(isfield(spin_system.control,{'catalog','ens_sizes','invariants','frozen_fields','worker_cases','pool_id'}))
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

