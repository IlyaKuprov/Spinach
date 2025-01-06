% A parallel wrapper around GRAPE that enables ensemble optimal control
% optimisations. This function handles systems with multiple control po-
% wer levels, multiple resonance offsets, multistate transfers, and en-
% sembles of drift Liouvillians. Syntax:
%
%      [fidelity,gradient,hessian]=ensemble(waveform,spin_system)
%
% Parameters:
%
%   waveform  - control coefficients for each control operator, rad/s
%
% Outputs:
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
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ensemble.m>

function [traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)

% Check consistency
grumble(spin_system,waveform);

% Get offset ensemble size
off_ens_sizes=cellfun(@numel,spin_system.control.offsets);
if ~isempty(off_ens_sizes)
    n_offset_vals=prod(off_ens_sizes);
else
    n_offset_vals=1;
end

% Extract ensemble grid dimensions
n_state_pairs=numel(spin_system.control.rho_init);     % State-target pair count
n_ens_systems=spin_system.control.ndrifts;             % Drift ensemble size
n_power_levls=numel(spin_system.control.pwr_levels);   % Power level count
n_phase_specs=size(spin_system.control.phase_cycle,1); % Phase cycle line count
n_cpm_mdepths=spin_system.control.cpm_mdp(3);          % Control power modulation depth count
n_cpm_phasept=spin_system.control.cpm_nph;             % Control power modulation phase count
n_distortions=size(spin_system.control.distortion,1);  % Distortion function ensemble size

% Create a catalog of the ensemble
catalog=(1:n_state_pairs)';
catalog=[kron(ones(n_ens_systems,1),catalog) kron((1:n_ens_systems)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_power_levls,1),catalog) kron((1:n_power_levls)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_offset_vals,1),catalog) kron((1:n_offset_vals)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_phase_specs,1),catalog) kron((1:n_phase_specs)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_cpm_mdepths,1),catalog) kron((1:n_cpm_mdepths)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_cpm_phasept,1),catalog) kron((1:n_cpm_phasept)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_distortions,1),catalog) kron((1:n_distortions)',ones(size(catalog,1),1))];

% Ensemble correlation: own state pair for each member
if ismember('rho_ens',spin_system.control.ens_corrs)
    catalog=catalog(:,2:end); 
    catalog=unique(catalog,'rows');
    catalog=[(1:size(catalog,1))' catalog];
end

% Ensemble correlation: own state pair for each drift
if ismember('rho_drift',spin_system.control.ens_corrs)
    catalog(catalog(:,1)~=catalog(:,2),:)=[];
end

% Ensemble correlation: own control power for each drift
if ismember('power_drift',spin_system.control.ens_corrs)
    catalog(catalog(:,3)~=catalog(:,2),:)=[];
end

% Count the cases
n_cases=size(catalog,1);

% Count the outputs
n_outputs=nargout;

% Preallocate outputs
traj_data=cell(n_cases,1); fidelities=cell(1,n_cases);
gradients=cell(1,n_cases); hessians=cell(1,n_cases);

% Waveform dimension statistics
ncont=size(waveform,1); nsteps=size(waveform,2);

% Parallel strategy
if strcmp(spin_system.control.parallel,'ensemble')

    % Over ensemble
    nworkers=poolsize;

elseif strcmp(spin_system.control.parallel,'time')

    % Pass down
    nworkers=0;

else

    % Complain and bomb out
    error('unknown parallelisation strategy.');

end
    
% Run the ensemble loop
parfor (n=1:n_cases,nworkers) %#ok<*PFBNS>
    
    % Extract ensemble indices
    n_rho=catalog(n,1); n_sys=catalog(n,2);
    n_pwr=catalog(n,3); n_off=catalog(n,4);
    n_phi=catalog(n,5); n_mdp=catalog(n,6);
    n_mph=catalog(n,7); n_dis=catalog(n,8);
    
    % Get initial and target state
    rho_init=spin_system.control.rho_init{n_rho};
    rho_targ=spin_system.control.rho_targ{n_rho};
    
    % Grab a copy of the waveform
    local_waveform=waveform;
    
    % Control power modulation depths
    mdepths=linspace(spin_system.control.cpm_mdp(1),...
                     spin_system.control.cpm_mdp(2),...
                     spin_system.control.cpm_mdp(3));

    % Current modulation depth
    mdepth=mdepths(n_mdp);

    % Control power modulation phases
    mphases=linspace(0,2*pi,spin_system.control.cpm_nph+1); 
    mphases=mphases(1:(end-1)); 
    
    % Current modulation phase
    mphase=mphases(n_mph);

    % Control power modulation function
    switch spin_system.control.integrator
        
        % Piecewise-constant
        case 'rectangle'

            % Equal number of points and intervals
            time_grid=cumsum(spin_system.control.pulse_dt)-...
                      spin_system.control.pulse_dt(1);
            amp_mod=1+mdepth*cos(2*pi*spin_system.control.cpm_frq*time_grid+mphase);
        
        % Piecewise-linear
        case 'trapezium'

            % One more point than intervals
            time_grid=[0 cumsum(spin_system.control.pulse_dt)];
            amp_mod=1+mdepth*cos(2*pi*spin_system.control.cpm_frq*time_grid+mphase);

        % Parfor needs this
        otherwise

            % Complain and bomb out
            error('unknown integrator.'); amp_mod=[]; %#ok<UNRCH> 

    end

    % Apply control power modulation
    local_waveform=local_waveform.*amp_mod;
    
    % Apply the phase cycle
    if ~isempty(spin_system.control.phase_cycle)
        
        % Apply phase to the initial state
        phi=spin_system.control.phase_cycle(n_phi,1);
        rho_init=exp(1i*phi)*rho_init;
        
        % Apply phase to the target state
        phi=spin_system.control.phase_cycle(n_phi,end);
        rho_targ=exp(1i*phi)*rho_targ;
        
        % Apply phases to the waveform
        for k=1:(size(local_waveform,1)/2)
            
            % Assemble complex waveform
            cplx_wave=local_waveform(2*k-1,:)+...
                   1i*local_waveform(2*k,:);
            
            % Get the phase
            phi=spin_system.control.phase_cycle(n_phi,k+1);
            
            % Apply the phase
            cplx_wave=exp(1i*phi)*cplx_wave;
            
            % Get back X and Y components
            local_waveform(2*k-1,:)=real(cplx_wave);
            local_waveform(2*k,:)=imag(cplx_wave);
            
        end
        
    end
    
    % Get drifts from pool ValueStore
    if (nworkers==0)||(~isworkernode)
        store=gcp('nocreate').ValueStore; 
        L=store(['oc_drift_' num2str(n_sys)]);
    else
        store=getCurrentValueStore(); 
        L=store(['oc_drift_' num2str(n_sys)]);
    end
      
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
            L=L+sparse(2*pi*spin_system.control.offsets{oper_idx}(vj)*...
                            spin_system.control.off_ops{oper_idx});

            % Next channel
            lin_idx=vi;

        end

    end

    % Move the waveform into physical units
    power_lvl=spin_system.control.pwr_levels(n_pwr);
    local_waveform=power_lvl*local_waveform;

    % Call GRAPE
    if n_outputs==2

        % Apply waveform distortions
        for k=1:size(spin_system.control.distortion,2)

            % Get distortion function
            dist_function=spin_system.control.distortion{n_dis,k};

            % Apply distortion function
            local_waveform=dist_function(local_waveform);

        end
            
        % Fidelity and trajectory
        switch spin_system.bas.formalism

            case {'sphten-liouv','zeeman-liouv'}
                
                % Call Liouville space version of the GRAPE function
                [traj_data{n},fidelities{n}]=grape_liouv(spin_system,L,spin_system.control.operators,...
                                                         local_waveform,rho_init,rho_targ,...
                                                         spin_system.control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj_data{n},fidelities{n}]=grape_hilb(spin_system,L,spin_system.control.operators,...
                                                        local_waveform,rho_init,rho_targ,...
                                                        spin_system.control.fidelity);
            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end
                                       
    elseif n_outputs==3

        % Get the Jacobian going
        J=speye(numel(local_waveform));

        % Apply waveform distortions
        for k=1:size(spin_system.control.distortion,2)

            % Get distortion function
            dist_function=spin_system.control.distortion{n_dis,k};

            % Apply distortion and get its Jacobian
            [local_waveform,stage_jacobian]=dist_function(local_waveform);

            % Combine Jacobians
            J=stage_jacobian*J;

        end
            
        % Fidelity and trajectory
        switch spin_system.bas.formalism

            case {'sphten-liouv','zeeman-liouv'}
                
                % Call Liouville space version of the GRAPE function
                [traj_data{n},fidelities{n},gradients{n}]=grape_liouv(spin_system,L,spin_system.control.operators,...
                                                                      local_waveform,rho_init,rho_targ,...
                                                                      spin_system.control.fidelity);
            case 'zeeman-hilb'

                % Call Hilbert space version of the GRAPE function
                [traj_data{n},fidelities{n},gradients{n}]=grape_hilb(spin_system,L,spin_system.control.operators,...
                                                                     local_waveform,rho_init,rho_targ,...
                                                                     spin_system.control.fidelity);
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
        switch spin_system.bas.formalism

            case {'sphten-liouv','zeeman-liouv'}
                
                % Call Liouville space version of the GRAPE function
                [traj_data{n},fidelities{n},...
                 gradients{n},hessians{n}]=grape_liouv(spin_system,L,spin_system.control.operators,...
                                                       local_waveform,rho_init,rho_targ,...
                                                       spin_system.control.fidelity);
            case 'zeeman-hilb'

                % Complain and bomb out
                error('Newton-Raphson methods are not available in Hilbert space, use LBFGS.');

            otherwise

                % Complain and bomb out
                error('unrecognised formalism specification.');

        end
                                                         
    end
    
    % Post-process gradient
    if (~isempty(spin_system.control.phase_cycle))&&(n_outputs>2)
        
        % Un-apply phases to gradient
        for k=1:(size(gradients{n},1)/2)
            
            % Assemble complex gradient
            cplx_grad=gradients{n}(2*k-1,:)+...
                   1i*gradients{n}(2*k,:);
            
            % Get the phase
            phi=spin_system.control.phase_cycle(n_phi,k+1);
            
            % Un-apply the phase
            cplx_grad=exp(-1i*phi)*cplx_grad;
            
            % Get back X and Y components
            gradients{n}(2*k-1,:)=real(cplx_grad);
            gradients{n}(2*k,:)=imag(cplx_grad);
            
        end
        
    end
    
    % Post-process Hessian
    if (~isempty(spin_system.control.phase_cycle))&&(n_outputs>3)
        
        % Re-shape the Hessian as [ncont x nsteps x nsteps x ncont]
        hessians{n}=reshape(hessians{n},[ncont nsteps ncont nsteps]);
    
        % Un-apply phases to Hessian
        for k=1:(size(gradients{n},1)/2)
            
            % Get the phase
            phi=spin_system.control.phase_cycle(n_phi,k+1);
            
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
    
    % Apply control power modulation
    if n_outputs>2
        gradients{n}=gradients{n}.*amp_mod;
    end
    if n_outputs>3
        hessians{n}=reshape(hessians{n},[ncont nsteps ncont nsteps]);
        hessians{n}=permute(hessians{n},[2 1 3 4]).*amp_mod';
        hessians{n}=permute(hessians{n},[4 1 3 2]).*amp_mod';
        hessians{n}=permute(hessians{n},[4 2 3 1]);
        hessians{n}=reshape(hessians{n},[ncont*nsteps nsteps*ncont]);
    else 
        hessians{n}=[];
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

% Run diagnostic plotting
if ~isempty(spin_system.control.plotting)
    
    % Send waveform, trajectories, and fidelities for plotting
    ctrl_trajan(spin_system,waveform,traj_data,fidelities);
    
end

end

% Consistency enforcement
function grumble(spin_system,waveform)
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be an array of real numbers.');
end
if size(waveform,1)~=numel(spin_system.control.operators)
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
        if strcmp(spin_system.bas.formalism,'zeeman-hilb')
            error('trapezium integration is not available in Hilbert space.');
        end
    otherwise
        error('unknown time propagation algorithm.');
end
end

% "After unsuccessful attempts to trap redtail monkeys
% at the Zika Forest with the intention of live-bleeding
% and release, monkeys had to be sampled by means of
% 12-bore shotguns."
%
% https://doi.org/10.1016/0035-9203(82)90161-4

