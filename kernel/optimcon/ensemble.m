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
%       own block, published by optimcon.m as pool constants; the per-
%       case physics runs in ens_block.m on each worker, so only the
%       waveform and the live control fields travel at each objective
%       evaluation, and the gradient, the Hessian, and averaged trajec-
%       tories are summed on the workers. This function must be called
%       from the client, on the pool that was open when optimcon.m ran:
%       a worker holds only its own block.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ensemble.m>

function [traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)

% Check consistency
grumble(spin_system,waveform,nargout);

% Worker-resident problem data handles
invariants=spin_system.control.invariants;
drift_slices=spin_system.control.drift_slices;

% Live problem data is the client-side control structure less what the workers already hold
control=rmfield(spin_system.control,{'invariants','drift_slices','worker_cases','basis'});
control.return_traj=isfield(control,'return_traj')&&control.return_traj;

% Count the outputs (fidelity always computed) and the cases
n_outputs=max(nargout,2); n_cases=size(control.catalog,1);

% Run the ensemble loop, each worker over its own case block
spmd (poolsize)

    % Evaluate the block of cases assigned to this worker
    [traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value,...
                                                           control,spmdIndex,waveform,n_outputs);

    % Reduce to the first worker and pack
    results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1),...
                   'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1));

end

% Collect from the first worker, fidelities back into catalog order
results=results{1}; gradient=results.grad; hessian=results.hess;
order=[spin_system.control.worker_cases{:}];
fidelities=zeros(1,n_cases); fidelities(order)=results.fid;

% Average the block trajectory sums, or put the trajectories back into catalog order
if ismember('average',control.traj_opts)
    ave_traj=results.traj{1}.forward;
    for n=2:numel(results.traj)
        ave_traj=ave_traj+results.traj{n}.forward;
    end
    traj_data={struct('forward',{(1/n_cases)*ave_traj})};
else
    traj_data=cell(n_cases,1); traj_data(order)=results.traj;
end

% Ensemble averages of fidelity, gradient, and Hessian
fidelity=sum(fidelities)/n_cases;
if n_outputs>2
    gradient=reshape(gradient/n_cases,size(waveform));
end
if n_outputs>3
    hessian=reshape(hessian/n_cases,numel(waveform)*[1 1]);
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

% Consistency enforcement
function grumble(spin_system,waveform,n_outputs)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if ~all(isfield(spin_system.control,{'catalog','ens_sizes','invariants','drift_slices','frozen_fields','worker_cases','pool_id'}))
    error('ensemble catalog missing from spin_system, run optimcon() first.');
end
current_pool=gcp('nocreate'); pool_id=0;
if ~isempty(current_pool), pool_id=current_pool.ID; end
if (~isempty(getCurrentWorker()))||(pool_id~=spin_system.control.pool_id)
    error('ensemble() must run on the client, on the pool that was open when optimcon() ran; re-run optimcon() after changing the pool.');
end
if any(isfield(spin_system.control,spin_system.control.frozen_fields))
    error('generators and operators are frozen after optimcon(), re-run optimcon() to change them.');
end
if (n_outputs>3)&&(~all(cellfun(@(f)isequal(f,@no_dist),spin_system.control.distortion(:))))
    error('Hessians are not available with waveform distortions.');
end
if (n_outputs>3)&&(~strcmp(spin_system.control.integrator,'rectangle'))
    error('Hessians are only available with the rectangle integrator.');
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

