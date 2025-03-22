% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian 
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. Uses Liouville-space formalism. Syntax:
%
%  [traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...
%                                       waveform,rho_init,rho_targ,...
%                                       fidelity_type)
% Parameters:
%
%   spin_system         - Spinach data object that has been through 
%                         the optimcon.m problem setup function.
% 
%   drifts              - the drift Liouvillians: a cell array con-
%                         taining one matrix (for time-independent 
%                         drift) or multiple matrices (for time-de-
%                         pendent drift).
%
%   controls            - control operators in Liouville space (cell 
%                         array of matrices).
%
%   waveform            - control coefficients for each control ope-
%                         rator (in vertical dimension) at each time
%                         step (in horizonal dimension), rad/s
%
%   rho_init            - initial state of the system as a vector in
%                         Liouville space.
%
%   rho_targ            - target state of the system as a vector in
%                         Liouville space.
%
%   fidelity_type       - 'real'   (real part of the overlap)
%                         'imag'   (imaginary part of the overlap)
%                         'square' (absolute square of the overlap)
%
% Outputs:
%
%   fidelity            - fidelity of the control sequence
%
%   grad                - gradient of the fidelity with respect to
%                         the control sequence
%
%   hess                - Hessian of the fidelity with respect to the 
%                         control sequence
%
%   traj_data.forward   - forward trajectory from the initial condi-
%                         tion(a stack of state vectors)
%
% Note: this is a low level function that is not designed to be called 
%       directly. Use grape_xy.m and grape_phase.m instead.
%
% david.goodwin@inano.au.dk
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% TODO (Keitel): add logic to avoid computing backward trajectory 
%                when the gradient is not requested
%
% <https://spindynamics.org/wiki/index.php?title=grape_liouv.m>

function [traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...
                                                    waveform,rho_init,rho_targ,...
                                                    fidelity_type) %#ok<*PFBNS>
% Check consistency
grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ);
    
% Count the outputs
n_outputs=nargout();

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Run array preallocations
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Number of time intervals and control operators
        nsteps=size(waveform,2); nctrls=size(waveform,1);

        % Preallocate forward and backward trajectories
        fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
        bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

        % Preallocate arrays used in Hessian calculation
        if n_outputs>3

            % Both Newton and Goodwin
            fwd_dP=cell(nctrls,nsteps);
            bwd_dP=cell(nctrls,nsteps); 
            fwd_d2P=cell(nctrls,nctrls,nsteps);

            % Goodwin method precomputes propagators
            if strcmp(spin_system.control.method,'goodwin')
                P=cell(1,nsteps);
            end

        else

            % Parfor needs these empty declarations
            fwd_dP={}; bwd_dP={}; fwd_d2P={}; P={};

        end

    % Piecewise-linear
    case 'trapezium'
      
        % Number of time intervals and control operators
        nsteps=size(waveform,2)-1; nctrls=size(waveform,1);

        % Preallocate forward and backward trajectories
        fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
        bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

        % Cannot do trapezium Hessians yet
        fwd_dP={}; bwd_dP={}; fwd_d2P={}; P={};
        
    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Hush up the output
spin_system.sys.output='hush';

% Pull the target back through the dead time using
% the last drift generator in the drift array
if spin_system.control.dead_time~=0
    rho_targ=step(spin_system,drifts{end}',rho_targ,...
                 -spin_system.control.dead_time);
end

% Push the source through the prefix sequence
% using the first element of the drift array
if ~isempty(spin_system.control.prefix)
    prefix=spin_system.control.prefix;
    rho_init=prefix(spin_system,drifts{1},rho_init);
end

% Push the target through the suffix sequence
% (which user needs to code in reverse time)
% using the last element of the drift array
if ~isempty(spin_system.control.suffix)
    suffix=spin_system.control.suffix;
    rho_targ=suffix(spin_system,drifts{end},rho_targ);
end

% Initialise forward and backward trajectories
fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ;

% Parallelisation prep
if n_outputs>2

    % Strip the spin system object for communication efficiency
    ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
    ss_parfor.bas.formalism=spin_system.bas.formalism;

    % Define a vector and a matrix of zeroes for auxiliary systems
    zero_state=complex(spalloc(size(rho_init,1),size(rho_init,2),0));
    zero_drift=complex(spalloc(size(drifts{1},1),size(drifts{1},2),0));

end

% Run forward and backward propagation
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Preallocate evolution generators
        L_forw=cell(1,nsteps); L_back=cell(1,nsteps);
        
        % Precompute evolution generators
        parfor n=1:nsteps

            % Decide current drifts
            if isscalar(drifts)

                % Time-independent drifts, including
                % conj-transpose for dissipative cases
                L_forw{n}=drifts{1}; L_back{n}=drifts{1}';

            else

                % Time-dependent drifts, including
                % conj-transpose for dissipative cases
                L_forw{n}=drifts{n}; L_back{n}=drifts{nsteps+1-n}';

            end

            % Add current controls to current drifts, including
            % conjugate-transpose for dissipative controls; the
            % waveform is always real
            for k=1:nctrls

                % Forward evolution generator
                L_forw{n}=L_forw{n}+waveform(k,n)*controls{k};

                % Backward evolution generator
                L_back{n}=L_back{n}+waveform(k,nsteps+1-n)*controls{k}';

            end

        end

        % Loop over time steps
        for n=1:nsteps

            % Get keyhole operators (projectors are Hermitian)
            keyhole_forw=spin_system.control.keyholes{n};
            keyhole_back=spin_system.control.keyholes{nsteps+1-n};

            % Apply keyhole to the backward trajectory
            if ~isempty(keyhole_back), bwd_traj(:,n)=keyhole_back(bwd_traj(:,n)); end

            % Goodwin's Hessian method pre-computes cumulative propagators
            if strcmp(spin_system.control.method,'goodwin')&&(n_outputs>3)
                P{n}=propagator(ss_parfor,L_forw{n},dt(n));
                if n>1
                    P{n}=P{n}*P{n-1};
                    P{n}=clean_up(spin_system,P{n},spin_system.tols.prop_chop);
                end
            end

            % Take a time step forwards and backwards
            fwd_traj(:,n+1)=step(spin_system,L_forw{n},fwd_traj(:,n),+dt(n));
            bwd_traj(:,n+1)=step(spin_system,L_back{n},bwd_traj(:,n),-dt(nsteps+1-n));

            % Apply keyhole to the forward trajectory
            if ~isempty(keyhole_forw), fwd_traj(:,n+1)=keyhole_forw(fwd_traj(:,n+1)); end

        end

    % Piecewise-linear
    case 'trapezium'

        % Preallocate evolution generators
        L_forw_left=cell(1,nsteps); L_forw_right=cell(1,nsteps);
        L_back_left=cell(1,nsteps); L_back_right=cell(1,nsteps);

        % Precompute evolution generators
        parfor n=1:nsteps

            % Decide current drifts
            if isscalar(drifts)

                % Time-independent drifts, including
                % conjugate-transpose for dissipative cases
                L_forw_left{n}=drifts{1};  L_forw_right{n}=drifts{1};
                L_back_left{n}=drifts{1}'; L_back_right{n}=drifts{1}';

            else

                % Time-dependent drifts, including
                % conjugate-transpose for dissipative cases
                L_forw_left{n}=drifts{n};           L_forw_right{n}=drifts{n+1};
                L_back_left{n}=drifts{nsteps+1-n}'; L_back_right{n}=drifts{nsteps+2-n}';

            end

            % Add current controls to current drifts, including
            % conjugate-transpose for dissipative controls; the
            % waveform is always real
            for k=1:nctrls
                L_forw_left{n}=L_forw_left{n}+waveform(k,n)*controls{k};
                L_forw_right{n}=L_forw_right{n}+waveform(k,n+1)*controls{k};
                L_back_left{n}=L_back_left{n}+waveform(k,nsteps+1-n)*controls{k}';
                L_back_right{n}=L_back_right{n}+waveform(k,nsteps+2-n)*controls{k}';
            end

        end

        % Loop over time steps
        for n=1:nsteps

            % Get keyhole operators (projectors are Hermitian)
            keyhole_forw=spin_system.control.keyholes{n};
            keyhole_back=spin_system.control.keyholes{nsteps+1-n};

            % Apply keyhole to the backward trajectory
            if ~isempty(keyhole_back), bwd_traj(:,n)=keyhole_back(bwd_traj(:,n)); end

            % Take a time step forwards and backwards
            fwd_traj(:,n+1)=step(spin_system,{ L_forw_left{n},...
                                              (L_forw_left{n}+L_forw_right{n})/2,...
                                               L_forw_right{n}},fwd_traj(:,n),+dt(n));
            bwd_traj(:,n+1)=step(spin_system,{ L_back_right{n},...
                                              (L_back_right{n}+L_back_left{n})/2,...
                                               L_back_left{n}},bwd_traj(:,n),-dt(nsteps+1-n));

            % Apply keyhole to the forward trajectory
            if ~isempty(keyhole_forw), fwd_traj(:,n+1)=keyhole_forw(fwd_traj(:,n+1)); end
            
        end

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Calculate the state overlap
overlap=rho_targ'*fwd_traj(:,end);

% Compute gradient
if n_outputs>2

    % Flip the backward trajectory
    bwd_traj=fliplr(bwd_traj);
    
    % Preallocate results
    grad=zeros(size(waveform),'like',1i);
    
    % Integrator-specific paths
    switch spin_system.control.integrator

        % Piecewise-constant
        case 'rectangle'

            % Loop over control sequence
            parfor n=1:nsteps
                
                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % Calculate gradient at this timestep
                for k=1:nctrls
            
                    % Create auxiliary system
                    aux_matrix=[ L_forw{n}   controls{k}
                                 zero_drift  L_forw{n}   ];

                    % Build the auxiliary vector
                    aux_vec=[zero_state; fwd_traj(:,n)];
            
                    % Propagate the auxiliary vector
                    aux_vec=step(ss_parfor,aux_matrix,aux_vec,dt(n));
            
                    % Compute the derivative
                    grad_col(k)=bwd_traj(:,n+1)'*aux_vec(1:(end/2));
            
                end
        
                % Add to the gradient array
                grad(:,n)=grad_col;

            end

        % Piecewise-linear
        case 'trapezium'

            % Pull the control commutators
            cc_comm_idx=spin_system.control.cc_comm_idx;
            cc_comm=spin_system.control.cc_comm;

            % Loop over control sequence
            parfor n=1:(nsteps+1)

                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % First step is special
                if n==1

                    % Decide current drift
                    if isscalar(drifts)

                        % Time-independent drift
                        L={drifts{1},drifts{1}};

                    else

                        % Time-dependent drift
                        L={drifts{1},drifts{2}};

                    end
    
                    % Loop over controls
                    for k=1:nctrls

                        % Build the auxiliary matrix
                        [DL_first,~]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n),waveform(:,1),waveform(:,2),k);

                        % Build the auxiliary vector
                        aux_vec=[zero_state; fwd_traj(:,n)];

                        % Propagate the auxiliary vector
                        aux_vec=step(ss_parfor,DL_first,aux_vec,dt(n));

                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,n+1)'*aux_vec(1:(end/2));

                    end
                
                % Last step is special
                elseif n==(nsteps+1)

                    % Decide current drift
                    if isscalar(drifts)

                        % Time-independent drift
                        L={drifts{1},drifts{1}};

                    else

                        % Time-dependent drift
                        L={drifts{n-1},drifts{n}};

                    end

                    % Loop over controls
                    for k=1:nctrls
                        
                        % Build the auxiliary matrix
                        [~,DR_last]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n-1),waveform(:,(end-1)),waveform(:,end),k);
                        
                        % Build the auxiliary vector
                        aux_vec=[zero_state; fwd_traj(:,n-1)];

                        % Propagate the auxiliary vector
                        aux_vec=step(ss_parfor,DR_last,aux_vec,dt(n-1));

                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,end)'*aux_vec(1:(end/2));

                    end

                % Middle steps
                else

                    % Loop over controls
                    for k=1:nctrls

                        % Left pair of drifts
                        if isscalar(drifts)

                            % Time-independent drift
                            L={drifts{1},drifts{1}};

                        else

                            % Time-dependent drift
                            L={drifts{n},drifts{n+1}};

                        end
                        
                        % Auxiliary matrix
                        [Right_DL,~]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n),waveform(:,n),waveform(:,n+1),k);

                        % Build the auxiliary vector
                        aux_vec_a=[zero_state; fwd_traj(:,n)];

                        % Propagate and extract derivative action
                        aux_vec_a=step(ss_parfor,Right_DL,aux_vec_a,dt(n)); 

                        % Product rule: [dP2]*[P1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n+1)'*aux_vec_a(1:(end/2));
                        
                        % Right pair of drifts
                        if isscalar(drifts)

                            % Time-independent drift
                            L={drifts{1},drifts{1}};

                        else

                            % Time-dependent drift
                            L={drifts{n-1},drifts{n}};

                        end

                        % Auxiliary vector and matrix
                        [~,Left_DR]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n-1),waveform(:,n-1),waveform(:,n),k);

                        % Build the auxiliary vector
                        aux_vec_b=[zero_state; fwd_traj(:,n-1)];
                        
                        % Propagate and extract derivative action
                        aux_vec_b=step(ss_parfor,Left_DR,aux_vec_b,dt(n));

                        % Product rule: [P2]*[dP1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n)'*aux_vec_b(1:(end/2));

                    end

                end

                % Add to the gradient array
                grad(:,n)=grad_col;

            end

        otherwise

            % Complain and bomb out
            error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

    end

end

% Compute Hessian
if strcmp(spin_system.control.integrator,'rectangle')&&(n_outputs>3)

    % Flip the backward trajectory for ease of paralellisation
    bwd_traj=fliplr(bwd_traj);

    % Compute derivative trajectories
    parfor n=1:nsteps

        % Preallocate local arrays
        fwd_dP_col=cell(nctrls,1);
        bwd_dP_col=cell(nctrls,1);
        fwd_d2P_block=cell(nctrls,nctrls);

        % Loop over control pairs
        for k=1:nctrls
            for j=1:nctrls

                % Diagonal elements
                if k==j

                    % Create forward auxiliary matrix
                    aux_matrix=[L_forw{n}     controls{k}   zero_drift
                                zero_drift    L_forw{n}     controls{j}
                                zero_drift    zero_drift    L_forw{n}  ];

                    % Create forward auxiliary vector
                    aux_vec=[zero_state; zero_state; fwd_traj(:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,dt(n));

                    % Only store acton by dP on rho
                    fwd_dP_col{k}=aux_vec((end/3+1):(2*end/3));

                    % Only store action by d2P on rho
                    fwd_d2P_block{k,j}=2*aux_vec(1:(end/3));

                    % Create backward auxiliary matrix
                    aux_matrix=[ L_back{n}   controls{k}'
                                 zero_drift  L_back{n}   ];

                    % Create backward auxiliary vector
                    aux_vec=[zero_state; bwd_traj(:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,-dt(nsteps+1-n));

                    % Only store the action by dP on rho
                    bwd_dP_col{k}=aux_vec(1:(end/2));

                % Off-diagonal elements
                else

                    % Create the forward auxiliary matrix
                    aux_matrix=[L_forw{n}     controls{k}   zero_drift
                                zero_drift    L_forw{n}     controls{j}
                                zero_drift    zero_drift    L_forw{n}  ];

                    % Create forward auxiliary vector
                    aux_vec=[zero_state; zero_state; fwd_traj(:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,dt(n));

                    % Only store action by d2P on rho
                    fwd_d2P_block{k,j}=2*aux_vec(1:(end/3));

                end

            end
        end

        % Store derivative trajectories
        fwd_dP(:,n)=fwd_dP_col; 
        bwd_dP(:,n)=bwd_dP_col;
        fwd_d2P(:,:,n)=fwd_d2P_block;

    end

    % Preallocate Hessian matrix
    hess=zeros(nctrls,nsteps,nctrls,nsteps,'like',1i);

    % Off-diagonal Hessian elements
    switch spin_system.control.method

        % Goodwin's method
        case 'goodwin'

            % Flip the backwards 
            % derivative trajectory
            bwd_dP=fliplr(bwd_dP);

            % Loop over timesteps
            parfor n=1:nsteps

                % Loop over controls
                for k=1:nctrls

                    % Propagate forward derivatives
                    % to first time step
                    fwd_dP{k,n}=P{n}'*fwd_dP{k,n};

                    % From second step
                    if n>1

                        % Propagate backward derivatives
                        % to first time step
                        bwd_dP{k,n}=bwd_dP{k,n}'*P{n-1};

                    end

                end

            end

            % Flip the backward trajectory
            bwd_traj=fliplr(bwd_traj);

            % Loop over timesteps
            parfor n=1:nsteps

                % Preallocate local Hessian column
                hess_col=zeros(nctrls,nsteps,nctrls,1,'like',1i);

                % Outer control loop
                for k=1:numel(controls)

                    % From second step
                    if n>1

                        % Inner control loop
                        for j=1:numel(controls)

                            % Construct array of forward derivatives
                            array_fwd_dP=cat(2,fwd_dP{j,1:n-1});

                            % Multiply out current backward derivatives and
                            % array of all forward derivatives
                            hess_col(j,1:n-1,k)=bwd_dP{k,n}*array_fwd_dP;

                        end

                    end

                    % Inner control loop
                    for j=1:numel(controls)

                        % Calculate non-mixed derivatives
                        hess_col(j,n,k)=bwd_traj(:,n+1)'*fwd_d2P{k,j,n};

                    end

                end

                % Add to Hessian array
                hess(:,:,:,n)=hess_col;

            end

        % Newton method
        case 'newton'

            % Flip the backward trajectory and backward derivatives
            bwd_traj=fliplr(bwd_traj); bwd_dP=fliplr(bwd_dP);

            % Loop over timesteps
            parfor n=1:nsteps

                % Allocate local Hessian column
                hess_col=zeros(nctrls,nsteps,nctrls,1,'like',1i);

                % Outer control loop
                for k=1:nctrls

                    % Pull backward propagated left derivative
                    bwd_dPk=bwd_dP{k,n};

                    % From second step
                    if n>1

                        % Inner control loop
                        for j=1:nctrls

                            % Calculate mixed derivatives
                            hess_col(j,n-1,k)=bwd_dPk'*fwd_dP{j,n-1};

                        end

                        % Loop over remaining time slices
                        for m=n-2:-1:1

                            % Decide current drifts
                            if isscalar(drifts)
                                L_back=drifts{1}';
                            else
                                L_back=drifts{m+1}';
                            end

                            % Inner control loop
                            for j=1:numel(controls)

                                % Add current controls to current drifts
                                L_back=L_back+waveform(j,m+1)*controls{j}';

                            end

                            % Propagate left derivative backwards
                            bwd_dPk=(step(spin_system,L_back,bwd_dPk,-dt(m+1)));

                            % Inner control loop
                            for j=1:numel(controls)

                                % Calculate mixed derivatives
                                hess_col(j,m,k)=bwd_dPk'*fwd_dP{j,m};

                            end

                        end

                    end

                    % Inner control loop
                    for j=1:numel(controls)

                        % Calculate non-mixed derivatives
                        hess_col(j,n,k)=bwd_traj(:,n+1)'*fwd_d2P{k,j,n};

                    end

                end

                % Add to Hessian array
                hess(:,:,:,n)=hess_col;

            end

        otherwise

            % Complain and bomb out
            error('Hessian calculation methods are ''newton'' and ''goodwin''')

    end

    % Merge the blocks and reorder
    hess=reshape(hess,nsteps*nctrls,nsteps*nctrls);
    
    % Force Hessian symmetry and fill empty Hessian entries
    hess=(hess+hess.').*(~kron(eye(nsteps),ones(nctrls,nctrls)))+...
         (hess+hess.').*( kron(eye(nsteps),ones(nctrls,nctrls)))/2;

end

% Fidelity and its derivatives
switch fidelity_type
    
    case {'real'}
        
        % Real part of the overlap
        fidelity=real(overlap);
        
        % Update Hessian
        if exist('hess','var'), hess=real(hess); end
        
        % Update gradient
        if exist('grad','var'), grad=real(grad); end
        
    case {'imag'}
        
        % Imaginary part of the overlap
        fidelity=imag(overlap);
        
        % Update Hessian
        if exist('hess','var'), hess=imag(hess); end
        
        % Update gradient
        if exist('grad','var'), grad=imag(grad); end
        
    case {'square'}
        
        % Absolute square of the overalp
        fidelity=overlap*conj(overlap);
        
        % Update Hessian
        if exist('hess','var')
            
            % Product rule
            hess=hess*conj(overlap)+grad(:)*transpose(conj(grad(:)))+...
                 conj(grad(:))*transpose(grad(:))+conj(hess)*overlap;
            
            % Cleaning up
            hess=real(hess);
        
        end
        
        % Update gradient
        if exist('grad','var')
            
            % Product rule
            grad=grad*conj(overlap)+overlap*conj(grad);
        
            % Cleaning up
            grad=real(grad);
        
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown fidelity type');
        
end

% Return trajectory data
traj_data.forward=fwd_traj;

% Catch unreachable objectives
if abs(fidelity)==0
    spin_system.sys.output=1;
    report(spin_system,'exactly zero fidelity: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end
if exist('grad','var')&&(norm(grad,1)==0)
    spin_system.sys.output=1;
    report(spin_system,'exactly zero gradient: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end

end

% Consistency enforcement
function grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv',...
                                        'zeeman-liouv',...
                                        'zeeman-wavef'})
    error('this function requires a state vector based formalism.');
end
if (~isnumeric(rho_init))||(~iscolumn(rho_init))
    error('rho_init must be a column vector.');
end
if (~isnumeric(rho_targ))||(~iscolumn(rho_init))
    error('rho_targ must be a column vector.');
end
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)     
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
    if (size(drifts{n},1)~=numel(rho_init))||...
       (size(drifts{n},1)~=numel(rho_targ))
        error('dimensions of drift, rho_init and rho_targ must be consistent.');
    end
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
       (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drifts{1},1))
        error('control operators must have the same size as drift operators.');
    end
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be a real numeric array.');
end
if size(waveform,1)~=numel(controls)
    error('number of waveform rows must be equal to the number of controls.');
end
if size(waveform,2)~=spin_system.control.pulse_ntpts
    error(['waveform must have ' int2str(spin_system.control.pulse_ntpts) ' columns.']);
end
if strcmp(spin_system.control.integrator,'rectangle')&&...
   (size(spin_system.control.pulse_dt,2)~=size(waveform,2))
    error('pulse_dt must have the same length as waveform for rectangle integrator');
end
if strcmp(spin_system.control.integrator,'trapezium')&&...
   (size(spin_system.control.pulse_dt,2)+1~=size(waveform,2))
    error('pulse_dt must be one element shorter than waveform for trapezium integrator');
end
end

% In any culture, subculture, or family in which belief is valued above
% thought, self-surrender is valued above self-expression, and conformity
% is valued above integrity, those who preserve their self-esteem are
% likely to be heroic exceptions.
%
% Nathaniel Branden

