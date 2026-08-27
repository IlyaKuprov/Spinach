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
%       are distributed over the parallel pool workers. Each case
%       fetches the frozen problem data from the pool constant pub-
%       lished by optimcon.m and grafts the live client-side control
%       structure on top of it, so only the waveform and the live
%       control fields travel at each objective evaluation.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ensemble.m>

function [traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)

% Check consistency
grumble(spin_system,waveform);

% Pull the ensemble case catalog
catalog=spin_system.control.catalog;
n_cases=size(catalog,1);

% Pull the worker-resident problem data handle
invariants=spin_system.control.invariants;

% Live problem data is the client-side control structure
control=rmfield(spin_system.control,'invariants');

% Default the trajectory return flag
control.return_traj=isfield(control,'return_traj')&&control.return_traj;

% Get offset ensemble size
off_ens_sizes=cellfun(@numel,control.offsets);

% Count the outputs
n_outputs=nargout;

% Preallocate outputs
traj_data=cell(n_cases,1); fidelities=cell(1,n_cases);
gradients=cell(1,n_cases); hessians=cell(1,n_cases);

% Waveform dimension statistics
ncont=size(waveform,1); nsteps=size(waveform,2);

% Parallelise over the ensemble
nworkers=poolsize;

% Run the ensemble loop
parfor (n=1:n_cases,nworkers) %#ok<*PFBNS>

    % Fetch worker-resident problem data
    ss=invariants.Value;

    % Graft live client data over the frozen worker copy
    frozen=ss.control; ss.control=control;
    for k=1:numel(control.frozen_fields)
        fname=control.frozen_fields{k};
        ss.control.(fname)=frozen.(fname);
    end

    % Extract ensemble indices
    n_rho=catalog(n,1); n_sys=catalog(n,2);
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
                [traj_data{n},fidelities{n}]=grape_liouv(ss,L,ss.control.operators,...
                                                         local_waveform,rho_init,rho_targ,...
                                                         control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj_data{n},fidelities{n}]=grape_hilb(ss,L,ss.control.operators,...
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
                [traj_data{n},fidelities{n},gradients{n}]=grape_liouv(ss,L,ss.control.operators,...
                                                                      local_waveform,rho_init,rho_targ,...
                                                                      control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj_data{n},fidelities{n},gradients{n}]=grape_hilb(ss,L,ss.control.operators,...
                                                                     local_waveform,rho_init,rho_targ,...
                                                                     control.fidelity);

            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end

        % Store the gradient layout
        [n_rows,n_cols]=size(gradients{n});

        % Stretch and apply the Jacobian
        gradients{n}=J'*gradients{n}(:);

        % Restore the original gradient layout
        gradients{n}=reshape(gradients{n},[n_rows n_cols]);

    elseif n_outputs==4

        % Fidelity and trajectory
        switch ss.bas.formalism

            case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}

                % Call Liouville space version of the GRAPE function
                [traj_data{n},fidelities{n},...
                 gradients{n},hessians{n}]=grape_liouv(ss,L,ss.control.operators,...
                                                       local_waveform,rho_init,rho_targ,...
                                                       control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj_data{n},fidelities{n},...
                 gradients{n},hessians{n}]=grape_hilb(ss,L,ss.control.operators,...
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
        for k=1:(size(gradients{n},1)/2)

            % Assemble complex gradient
            cplx_grad=gradients{n}(2*k-1,:)+...
                   1i*gradients{n}(2*k,:);

            % Get the phase
            phi=control.phase_cycle(n_phi,k+1);

            % Un-apply the phase
            cplx_grad=exp(-1i*phi)*cplx_grad;

            % Get back X and Y components
            gradients{n}(2*k-1,:)=real(cplx_grad);
            gradients{n}(2*k,:)=imag(cplx_grad);

        end

    end

    % Post-process Hessian
    if (~isempty(control.phase_cycle))&&(n_outputs>3)

        % Re-shape the Hessian as [ncont x nsteps x nsteps x ncont]
        hessians{n}=reshape(hessians{n},[ncont nsteps ncont nsteps]);

        % Un-apply phases to Hessian
        for k=1:(size(gradients{n},1)/2)

            % Get the phase
            phi=control.phase_cycle(n_phi,k+1);

            % Assemble complex Hessian - left
            cplx_hess=hessians{n}(2*k-1,:,:,:)+...
                   1i*hessians{n}(2*k,:,:,:);

            % Un-apply the phase
            cplx_hess=exp(-1i*phi)*cplx_hess;

            % Get back X and Y components
            hessians{n}(2*k-1,:,:,:)=real(cplx_hess);
            hessians{n}(2*k,:,:,:)=imag(cplx_hess);

            % Assemble complex Hessian - right
            cplx_hess=hessians{n}(:,:,2*k-1,:)+...
                   1i*hessians{n}(:,:,2*k,:);

            % Un-apply the phase
            cplx_hess=exp(-1i*phi)*cplx_hess;

            % Get back X and Y components
            hessians{n}(:,:,2*k-1,:)=real(cplx_hess);
            hessians{n}(:,:,2*k,:)=imag(cplx_hess);

        end

        % Reshape the Hessian back
        hessians{n}=reshape(hessians{n},[ncont*nsteps nsteps*ncont]);

    end

    % Apply power level
    if n_outputs>2
        gradients{n}=power_lvl*gradients{n}(:);
    end
    if n_outputs>3
        hessians{n}=power_lvl*power_lvl*hessians{n}(:);
    end

end

% Apply trajectory options
if ismember('average',spin_system.control.traj_opts)

    % Return average trajectory
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

% Add up gradients
if n_outputs>2
    gradient=sum(cell2mat(gradients),2)/n_cases;
    gradient=reshape(gradient,size(waveform));
end

% Add up Hessians
if (n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')
    hessian=sum(cell2mat(hessians),2)/n_cases;
    hessian=reshape(hessian,numel(waveform)*[1 1]);
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

% Consistency enforcement
function grumble(spin_system,waveform)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if ~all(isfield(spin_system.control,{'catalog','ens_sizes','invariants','frozen_fields'}))
    error('ensemble catalog missing from spin_system, run optimcon() first.');
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

