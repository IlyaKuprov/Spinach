% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian 
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. Uses Liouville-space formalism. Syntax:
%
%        [traj_data,fidelity,...
%         grad,hess]=grape_liouv(spin_system,drifts,controls,...
%                                waveform,rho_init,rho_targ,...
%                                fidelity_type)
% Parameters:
%
%   spin_system         - Spinach data object that has been through 
%                         the optimcon.m problem setup function.
% 
%   drifts              - the drift Liouvillians: a cell array con-
%                         taining one matrix (for time-independent 
%                         drift) or multiple matrices (one per time
%                         slice / point, for time-dependent drift).
%
%   controls            - control operators in Liouville space (cell 
%                         array of matrices).
%
%   waveform            - control coefficients for each control ope-
%                         rator (in vertical dimension) at each time
%                         slice / point (horizonal dimension), rad/s
%
%   rho_init            - initial state of the system as a vector in
%                         Liouville space, ignored in stroboscopic
%                         steady state optimisations
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
%   hess                - Hessian of the fidelity with respect to 
%                         the control sequence, not available for
%                         piecewise-linear and stroboscopic stea-
%                         dy state optimisations
%
%   traj_data.forward   - forward trajectory from the initial con-
%                         dition or stroboscopic steady state (a 
%                         stack of state vectors)
%
% Note: this is a low level function that is not designed to be called 
%       directly. Use grape_xy.m, grape_phase.m, or other wrapper func-
%       tions instead.
%
% Note: trajectory cost terms are read from spin_system.control: when
%       fid_type is 'average', the fidelity is averaged over the pulse
%       nodes 1..N instead of being taken at the last node; traj_pen
%       operators are summed, their expectation value is averaged over
%       the same nodes and subtracted from the fidelity. Both terms use
%       costates that ride on the backward sweep, the trajectory never
%       leaves the worker. Hessians are not available with these terms.
%
% david.goodwin@inano.au.dk
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grape_liouv.m>

function [traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...
                                                    waveform,rho_init,rho_targ,...
                                                    fidelity_type)
% Check consistency
grumble(spin_system,drifts,controls,waveform,...
        rho_init,rho_targ,fidelity_type);
    
% Count the outputs
n_outputs=nargout();

% Pull the trajectory cost term settings
fid_avg=strcmp(spin_system.control.fid_type,'average');
pen_on=~isempty(spin_system.control.traj_pen);

% Trajectory cost terms have no Hessians
if (n_outputs>3)&&(fid_avg||pen_on)
    error('Hessians are not available with trajectory cost terms.');
end

% Trajectory cost terms need a waveform-independent initial state
if spin_system.control.steady&&(fid_avg||pen_on)
    error('trajectory cost terms are not available with stroboscopic steady states.');
end

% Phase cycle factors cancel in the overlap but not in the penalty
if pen_on&&(~isempty(spin_system.control.phase_cycle))
    error('trajectory penalties are not available with phase cycles.');
end

% Sum the trajectory penalty operators
if pen_on
    pen_op=spin_system.control.traj_pen{1};
    for k=2:numel(spin_system.control.traj_pen)
        pen_op=pen_op+spin_system.control.traj_pen{k};
    end
end

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Make freeze mask explicit
if isempty(spin_system.control.freeze)
    frozen=false(size(waveform));
else
    frozen=spin_system.control.freeze;
end

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

            % Goodwin Hessian route needs
            % cumulative propagators
            if strcmp(spin_system.control.method,'goodwin')
                P_cum=cell(1,nsteps);
            end

        else

            % Empty declarations tested downstream
            fwd_dP={}; bwd_dP={}; fwd_d2P={}; P_cum={};

        end

        % Steady state needs propagators
        if spin_system.control.steady
            P=cell(1,nsteps);
        else
            P={};
        end

    % Piecewise-linear
    case 'trapezium'
      
        % Number of time intervals and control operators
        nsteps=size(waveform,2)-1; nctrls=size(waveform,1);

        % Preallocate forward and backward trajectories
        fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
        bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

        % No Hessians or precomputed propagators for trapezium
        fwd_dP={}; bwd_dP={}; fwd_d2P={}; P_cum={}; P={};
        
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

% Define a vector and a matrix of zeroes for auxiliary systems
zero_state=complex(spalloc(size(rho_init,1),size(rho_init,2),0));
zero_drift=complex(spalloc(size(drifts{1},1),size(drifts{1},2),0));

% Initialise forward and backward trajectories
fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ;

% Count the drifts
ndrifts=numel(drifts);

% Pull Bloch-Siegert response operators
bss_on=isfield(spin_system.control,'bsiegert')&&spin_system.control.bsiegert;
if bss_on
    resp_ops=spin_system.control.resp_ops;
end

% Run forward and backward propagation
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Preallocate evolution generators
        L_forw=cell(1,nsteps); L_back=cell(1,nsteps);

        % Precompute evolution generators
        for n=1:nsteps

            % Cycle through the drifts array
            L_forw{n}=drifts{mod(n-1,ndrifts)+1};
            L_back{n}=drifts{mod(nsteps-n,ndrifts)+1}';

            % Add current controls to current drifts, including
            % conjugate-transpose for dissipative controls; the
            % waveform is always real
            for k=1:nctrls

                % Forward evolution generator
                L_forw{n}=L_forw{n}+waveform(k,n)*controls{k};

                % Backward evolution generator
                L_back{n}=L_back{n}+waveform(k,nsteps+1-n)*controls{k}';

            end

            % Add Bloch-Siegert terms to the generators
            if bss_on
                for k=1:nctrls
                    L_forw{n}=L_forw{n}+waveform(k,n)^2*resp_ops{k};
                    L_back{n}=L_back{n}+waveform(k,nsteps+1-n)^2*resp_ops{k}';
                end
            end

        end

        % Set the stage for StStSt
        if spin_system.control.steady

            % Compute all propagators
            P_tot=speye(size(rho_init,1));
            for n=1:nsteps

                % Get the step propagator
                P{n}=propagator(spin_system,L_forw{n},dt(n));

                % Merge into the total (physical: no keyholes)
                P_tot=clean_up(spin_system,P{n}*P_tot,...
                               spin_system.tols.prop_chop);

            end

            % Overwrite the initial state
            rho_init=steady(spin_system,P_tot);

            % Start trajectory
            fwd_traj(:,1)=rho_init;

            % Make sure the target state is a valid observable
            if rho_targ(1)~=0 
                error('target state must have a zero trace.'); 
            end

            % Tiptoe around the unit state singularity
            Q=(speye(size(P_tot))-P_tot)'; Q=Q(2:end,2:end);
            rho_targ_dressed=zeros(size(rho_targ),'like',1i);
            rho_targ_dressed(2:end)=Q\rho_targ(2:end);
            
            % Check if destination state dressing succeeded
            if any(~isfinite(rho_targ_dressed))
                error('destination state dressing failed.');
            end

            % Overwrite backward traj start
            bwd_traj(:,1)=rho_targ_dressed;
             
        end

    % Piecewise-linear
    case 'trapezium'

        % Preallocate evolution generators
        L_forw=cell(1,nsteps); L_back=cell(1,nsteps);

        % Precompute evolution generators
        for n=1:nsteps

            % Cycle through the drifts array
            forw_left=drifts{mod(n-1,ndrifts)+1};
            forw_right=drifts{mod(n,ndrifts)+1};
            back_left=drifts{mod(nsteps-n,ndrifts)+1}';
            back_right=drifts{mod(nsteps+1-n,ndrifts)+1}';

            % Add current controls to the node generators
            for k=1:nctrls
                forw_left=forw_left+waveform(k,n)*controls{k};
                forw_right=forw_right+waveform(k,n+1)*controls{k};
                back_left=back_left+waveform(k,nsteps+1-n)*controls{k}';
                back_right=back_right+waveform(k,nsteps+2-n)*controls{k}';
            end

            % Assemble the product quadrature generator triplets
            L_forw{n}={forw_left,(forw_left+forw_right)/2,forw_right};
            L_back{n}={back_right,(back_right+back_left)/2,back_left};

        end

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Run the forward trajectory
for n=1:nsteps

    % Get the forward keyhole operator
    keyhole_forw=spin_system.control.keyholes{n};

    % Goodwin's optimiser needs cumulative propagators
    if strcmp(spin_system.control.method,'goodwin')&&(n_outputs>3)
        P_cum{n}=propagator(spin_system,L_forw{n},dt(n));
        if n>1
            P_cum{n}=P_cum{n}*P_cum{n-1};
            P_cum{n}=clean_up(spin_system,P_cum{n},...
                              spin_system.tols.prop_chop);
        end
    end

    % Take a time step forwards
    if (~isempty(P))&&(~isempty(P{n}))
        fwd_traj(:,n+1)=P{n}*fwd_traj(:,n);
    else
        fwd_traj(:,n+1)=step(spin_system,L_forw{n},fwd_traj(:,n),+dt(n));
    end

    % Apply keyhole to forward trajectory
    if ~isempty(keyhole_forw)
        fwd_traj(:,n+1)=keyhole_forw(fwd_traj(:,n+1));
    end

end

% Overlaps with the target at the nodes
overlaps=rho_targ'*fwd_traj(:,2:end); overlap=overlaps(end);

% Trajectory penalty value and costate sources at the nodes
if pen_on
    if strcmp(spin_system.bas.formalism,'zeeman-wavef')
        pen_src=pen_op*fwd_traj(:,2:end); src_idx=1:nsteps;
        pen_val=mean(real(dot(pen_src,fwd_traj(:,2:end)))); pen_src=2*pen_src/nsteps;
    else
        pen_val=mean(real(pen_op'*fwd_traj(:,2:end)));
        pen_src=pen_op/nsteps; src_idx=ones(1,nsteps);
    end
end

% Run the backward trajectories
if n_outputs>2

    % Fidelity costate source weights at the nodes
    if fid_avg&&strcmp(fidelity_type,'square')
        weights=overlaps/nsteps;
    elseif fid_avg
        weights=ones(1,nsteps)/nsteps;
    else
        weights=[zeros(1,nsteps-1) 1];
    end

    % Scale the fidelity costate start by the last node weight
    bwd_traj(:,1)=weights(nsteps)*bwd_traj(:,1);

    % Start the penalty costate at the last node
    if pen_on
        pen_traj=zeros(size(fwd_traj),'like',1i);
        pen_traj(:,1)=pen_src(:,src_idx(nsteps));
    end

    % Loop over time steps
    for n=1:nsteps

        % Get the backward keyhole operator
        keyhole_back=spin_system.control.keyholes{nsteps+1-n};

        % Apply keyhole to the costates
        if ~isempty(keyhole_back)
            bwd_traj(:,n)=keyhole_back(bwd_traj(:,n));
            if pen_on, pen_traj(:,n)=keyhole_back(pen_traj(:,n)); end
        end

        % Take a time step backwards
        if (~isempty(P))&&(~isempty(P{n}))
            bwd_traj(:,n+1)=P{nsteps+1-n}'*bwd_traj(:,n);
        else
            bwd_traj(:,n+1)=step(spin_system,L_back{n},bwd_traj(:,n),-dt(nsteps+1-n));
        end
        if pen_on
            pen_traj(:,n+1)=step(spin_system,L_back{n},pen_traj(:,n),-dt(nsteps+1-n));
        end

        % Add the sources at the previous node
        if n<nsteps
            bwd_traj(:,n+1)=bwd_traj(:,n+1)+weights(nsteps-n)*rho_targ;
            if pen_on, pen_traj(:,n+1)=pen_traj(:,n+1)+pen_src(:,src_idx(nsteps-n)); end
        end

    end

    % Flip the costates to match forward indexing
    bwd_traj=fliplr(bwd_traj);
    if pen_on, pen_traj=fliplr(pen_traj); end

end

% Compute gradient
if n_outputs>2

    % Preallocate results
    grad=zeros(size(waveform),'like',1i);
    if pen_on, pen_grad=zeros(size(waveform)); end
    
    % Integrator-specific paths
    switch spin_system.control.integrator

        % Piecewise-constant
        case 'rectangle'

            % Over time steps
            for n=1:nsteps

                % Preallocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % Over channels
                for k=1:nctrls

                    % Check the freeze mask
                    if ~frozen(k,n)

                        % Pull the derivative direction through the Bloch-Siegert map
                        if bss_on
                            deriv_op=controls{k}+2*waveform(k,n)*resp_ops{k};
                        else
                            deriv_op=controls{k};
                        end

                        % Create auxiliary system
                        aux_matrix=[ L_forw{n}   deriv_op
                                     zero_drift  L_forw{n}   ];

                        % Build the auxiliary vector
                        aux_vec=[zero_state; fwd_traj(:,n)];
            
                        % Propagate the auxiliary vector
                        aux_vec=step(spin_system,aux_matrix,aux_vec,dt(n));
                
                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,n+1)'*aux_vec(1:(end/2));
                        if pen_on, pen_grad(k,n)=real(pen_traj(:,n+1)'*aux_vec(1:(end/2))); end

                    else

                        % No step
                        grad_col(k)=0;

                    end
            
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
            for n=1:(nsteps+1)

                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % First step is special
                if n==1

                    % Left pair of drifts
                    L={drifts{mod(n-1,ndrifts)+1},...
                       drifts{mod(n,ndrifts)+1}};
    
                    % Loop over controls
                    for k=1:nctrls

                        % Build the auxiliary matrix
                        [DL_first,~]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n),waveform(:,1),waveform(:,2),k);

                        % Build the auxiliary vector
                        aux_vec=[zero_state; fwd_traj(:,n)];

                        % Propagate the auxiliary vector
                        aux_vec=step(spin_system,DL_first,aux_vec,dt(n));

                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,n+1)'*aux_vec(1:(end/2));
                        if pen_on, pen_grad(k,n)=real(pen_traj(:,n+1)'*aux_vec(1:(end/2))); end

                    end
                
                % Last step is special
                elseif n==(nsteps+1)

                    % Right pair of drifts
                    L={drifts{mod(n-2,ndrifts)+1},...
                       drifts{mod(n-1,ndrifts)+1}};
                    
                    % Loop over controls
                    for k=1:nctrls
                        
                        % Build the auxiliary matrix
                        [~,DR_last]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n-1),waveform(:,(end-1)),waveform(:,end),k);
                        
                        % Build the auxiliary vector
                        aux_vec=[zero_state; fwd_traj(:,n-1)];

                        % Propagate the auxiliary vector
                        aux_vec=step(spin_system,DR_last,aux_vec,dt(n-1));

                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,end)'*aux_vec(1:(end/2));
                        if pen_on, pen_grad(k,n)=real(pen_traj(:,end)'*aux_vec(1:(end/2))); end

                    end

                % Middle steps
                else

                    % Loop over controls
                    for k=1:nctrls

                        % Left pair of drifts
                        L={drifts{mod(n-1,ndrifts)+1},...
                           drifts{mod(n,ndrifts)+1}};
                        
                        % Auxiliary matrix
                        [Right_DL,~]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n),waveform(:,n),waveform(:,n+1),k);

                        % Build the auxiliary vector
                        aux_vec_a=[zero_state; fwd_traj(:,n)];

                        % Propagate and extract derivative action
                        aux_vec_a=step(spin_system,Right_DL,aux_vec_a,dt(n)); 

                        % Product rule: [dP2]*[P1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n+1)'*aux_vec_a(1:(end/2));
                        if pen_on, pen_grad(k,n)=pen_grad(k,n)+real(pen_traj(:,n+1)'*aux_vec_a(1:(end/2))); end
                        
                        % Right pair of drifts
                        L={drifts{mod(n-2,ndrifts)+1},...
                           drifts{mod(n-1,ndrifts)+1}};

                        % Auxiliary vector and matrix
                        [~,Left_DR]=aux_mat(L,controls,cc_comm_idx,cc_comm,dt(n-1),waveform(:,n-1),waveform(:,n),k);

                        % Build the auxiliary vector
                        aux_vec_b=[zero_state; fwd_traj(:,n-1)];
                        
                        % Propagate and extract derivative action
                        aux_vec_b=step(spin_system,Left_DR,aux_vec_b,dt(n-1));

                        % Product rule: [P2]*[dP1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n)'*aux_vec_b(1:(end/2));
                        if pen_on, pen_grad(k,n)=pen_grad(k,n)+real(pen_traj(:,n)'*aux_vec_b(1:(end/2))); end

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

    % Flip the backward trajectory to match forward indexing
    bwd_traj=fliplr(bwd_traj);

    % Compute derivative trajectories
    for n=1:nsteps

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
            for n=1:nsteps

                % Loop over controls
                for k=1:nctrls

                    % Propagate forward derivatives
                    % to first time step
                    fwd_dP{k,n}=P_cum{n}'*fwd_dP{k,n};

                    % From second step
                    if n>1

                        % Propagate backward derivatives
                        % to first time step
                        bwd_dP{k,n}=bwd_dP{k,n}'*P_cum{n-1};

                    end

                end

            end

            % Flip the backward trajectory
            bwd_traj=fliplr(bwd_traj);

            % Loop over timesteps
            for n=1:nsteps

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
            for n=1:nsteps

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

                            % Loop through the drifts array
                            L_back=drifts{mod(m,ndrifts)+1}';

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
        
        % Real part of the overlap, averaged over the nodes if requested
        if fid_avg
            fidelity=mean(real(overlaps));
        else
            fidelity=real(overlap);
        end
        
        % Update Hessian
        if exist('hess','var'), hess=real(hess); end
        
        % Update gradient
        if exist('grad','var'), grad=real(grad); end
        
    case {'imag'}
        
        % Imaginary part of the overlap, averaged over the nodes if requested
        if fid_avg
            fidelity=mean(imag(overlaps));
        else
            fidelity=imag(overlap);
        end
        
        % Update Hessian
        if exist('hess','var'), hess=imag(hess); end
        
        % Update gradient
        if exist('grad','var'), grad=imag(grad); end
        
    case {'square'}
        
        % Absolute square of the overlap, averaged over the nodes if requested
        if fid_avg
            fidelity=mean(abs(overlaps).^2);
        else
            fidelity=overlap*conj(overlap);
        end
        
        % Update Hessian
        if exist('hess','var')
            
            % Product rule
            hess=hess*conj(overlap)+...
                 grad(:)*transpose(conj(grad(:)))+...
                 conj(grad(:))*transpose(grad(:))+...
                 conj(hess)*overlap;
            
            % Cleaning up
            hess=real(hess);

        end
        
        % Update gradient
        if exist('grad','var')
            
            % Product rule, averaged overlaps are inside the costate
            if fid_avg
                grad=2*real(grad);
            else
                grad=real(grad*conj(overlap)+overlap*conj(grad));
            end
        
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown fidelity type');
        
end

% Decouple frozen directions
if exist('hess','var')
    hess(frozen(:),:)=0; 
    hess(:,frozen(:))=0; 
    hess(frozen(:),frozen(:))=1;
end

% Return the trajectory (a huge array) only if needed
if (isfield(spin_system.control,'return_traj')&&spin_system.control.return_traj)||...
   any(ismember({'correlation_order','coherence_order',...
                 'local_each_spin',  'total_each_spin',...
                 'level_populations','trajectory'},spin_system.control.plotting(:)))
    traj_data.forward=fwd_traj;
else
    traj_data.forward=[];
end

% Catch unreachable terminal objectives, trajectory cost terms may cancel legitimately
if (~(fid_avg||pen_on))&&(abs(fidelity)==0)
    spin_system.sys.output=1;
    report(spin_system,'exactly zero fidelity: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end
if (~(fid_avg||pen_on))&&exist('grad','var')&&(norm(grad,1)==0)
    spin_system.sys.output=1;
    report(spin_system,'exactly zero gradient: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end

% Subtract the trajectory penalty
if pen_on
    fidelity=fidelity-pen_val;
    if exist('grad','var'), grad=grad-pen_grad; end
end

end

% Consistency enforcement
function grumble(spin_system,drifts,controls,waveform,...
                 rho_init,rho_targ,fidelity_type)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv',...
                                        'zeeman-liouv',...
                                        'zeeman-wavef'})
    error('this function requires a state vector based formalism.');
end
if isfield(spin_system.control,'steady')&&spin_system.control.steady
    if ismember(spin_system.control.method,{'newton','goodwin'})
        error('Newton-Raphson unavailable for stroboscopic steady states.');
    end
    if ~strcmp(spin_system.control.integrator,'rectangle')
        error('only rectangle integrator supported for stroboscopic steady states.');
    end
    if spin_system.control.dead_time~=0
        error('dead times are inapplicable to stroboscopic steady states.');
    end
    if ~isempty(spin_system.control.prefix) || ...
       ~isempty(spin_system.control.suffix)
        error('prefixes and suffixes are inapplicable to stroboscopic steady states.');
    end
    if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
        error('sphten-liouv formalism is required for stroboscopic steady states.');
    end
end
if (~isnumeric(rho_init))||(~iscolumn(rho_init))
    error('rho_init must be a column vector.');
end
if (~isnumeric(rho_targ))||(~iscolumn(rho_targ))
    error('rho_targ must be a column vector.');
end
if (~ischar(fidelity_type))||(~ismember(fidelity_type,{'real','imag','square'}))
    error('fidelity_type must be ''real'', ''imag'', or ''square''.');
end
if (~ischar(spin_system.control.fid_type))||...
   (~ismember(spin_system.control.fid_type,{'terminal','average'}))
    error('spin_system.control.fid_type must be ''terminal'' or ''average''.');
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

