% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. Uses Hilbert-space formalism. Syntax:
%
%  [traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...
%                                       waveform,rho_init,rho_targ,...
%                                       fidelity_type)
% Parameters:
%
%   spin_system         - Spinach data object that has been through
%                         the optimcon.m problem setup function.
%
%   drifts              - the drift Hamiltonians: a cell array con-
%                         taining one matrix (for time-independent
%                         drift) or multiple matrices (for time-de-
%                         pendent drift).
%
%   controls            - control operators in Hilbert space (cell
%                         array of matrices).
%
%   waveform            - control coefficients for each control ope-
%                         rator (in vertical dimension) at each time
%                         step (in horizonal dimension), rad/s
%
%   rho_init            - initial state of the system as a density
%                         matrix in Hilbert space.
%
%   rho_targ            - target state of the system as a density
%                         matrix in Hilbert space.
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
%                         tion(a stack of state matrices)
%
% Note: this is a low level function that is not designed to be called
%       directly. Use grape_xy.m and grape_phase.m instead.
%
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grape_hilb.m>

function [traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...
                                                   waveform,rho_init,rho_targ,...
                                                   fidelity_type) %#ok<*PFBNS>

% Check consistency
grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type);

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

    % Piecewise-linear
    case 'trapezium'

        % Number of time intervals and control operators
        nsteps=size(waveform,2)-1; nctrls=size(waveform,1);

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Hush up the output
ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
ss_parfor.bas.formalism=spin_system.bas.formalism;
ss_parfor.sys.output='hush';

% Pull the target back through the dead time using
% the last drift generator in the drift array
if spin_system.control.dead_time~=0
    rho_targ=step(spin_system,drifts{end},rho_targ,...
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

% Preallocate forward trajectory
fwd_traj=cell(1,nsteps+1); fwd_traj{1}=rho_init;

% Preallocate backward trajectory
if n_outputs>2
    bwd_traj=cell(1,nsteps+1); bwd_traj{1}=rho_targ;
end

% Count the drifts
ndrifts=numel(drifts);

% Precompute interval Hamiltonians and propagators
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Preallocate evolution generators and propagators
        H=cell(1,nsteps); P=cell(1,nsteps);

        % Precompute evolution generators and propagators
        parfor n=1:nsteps

            % Cycle through the drifts array
            H{n}=drifts{mod(n-1,ndrifts)+1};

            % Add current controls
            for k=1:nctrls

                % Forward evolution generator
                H{n}=H{n}+waveform(k,n)*controls{k};

            end

            % Compute the interval propagator
            P{n}=propagator(ss_parfor,H{n},dt(n));

        end

    % Piecewise-linear
    case 'trapezium'

        % Preallocate interval generators and propagators
        H=cell(1,nsteps); P=cell(1,nsteps);

        % Precompute edge generators and propagators
        parfor n=1:nsteps

            % Cycle through the drifts array
            left_ham=drifts{mod(n-1,ndrifts)+1};
            right_ham=drifts{mod(n,ndrifts)+1};

            % Add current controls to the interval edges
            for k=1:nctrls
                left_ham=left_ham+waveform(k,n)*controls{k};
                right_ham=right_ham+waveform(k,n+1)*controls{k};
            end

            % Build the Iserles-Norsett product-quadrature generator
            H{n}=(left_ham+right_ham)/2+...
                 1i*dt(n)*(1/12)*(left_ham*right_ham-...
                                   right_ham*left_ham);

            % Compute the interval propagator
            P{n}=propagator(ss_parfor,H{n},dt(n));

        end

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Run the forward trajectory
for n=1:nsteps

    % Take a time step forwards
    fwd_traj{n+1}=P{n}*fwd_traj{n}*P{n}';

    % Get the forward keyhole operator
    keyhole_forw=spin_system.control.keyholes{n};

    % Apply keyhole to the forward trajectory
    if ~isempty(keyhole_forw), fwd_traj{n+1}=keyhole_forw(fwd_traj{n+1}); end

end

% Run the backward trajectory
if n_outputs>2
    for n=1:nsteps

        % Index the backward interval
        m=nsteps+1-n;

        % Get the backward keyhole operator
        keyhole_back=spin_system.control.keyholes{m};

        % Apply keyhole to the backward trajectory
        if ~isempty(keyhole_back), bwd_traj{n}=keyhole_back(bwd_traj{n}); end

        % Take a time step backwards
        bwd_traj{n+1}=P{m}'*bwd_traj{n}*P{m};

    end

    % Flip the backward trajectory
    bwd_traj=fliplr(bwd_traj);

end

% Calculate the state overlap
overlap=hdot(fwd_traj{end},rho_targ);

% Compute gradient
if n_outputs>2

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

                    % Compute directional derivative of the propagator
                    auxmat=dirdiff(ss_parfor,H{n},controls{k},dt(n),2);

                    % Apply the derivative to the density matrix
                    rho_deriv=auxmat{2}*fwd_traj{n}*auxmat{1}'+...
                              auxmat{1}*fwd_traj{n}*auxmat{2}';

                    % Compute the derivative of the objective
                    grad_col(k)=hdot(bwd_traj{n+1},rho_deriv);

                end

                % Add to the gradient array
                grad(:,n)=grad_col;

            end

        % Piecewise-linear
        case 'trapezium'

            % Pull the control commutators
            cc_comm_idx=spin_system.control.cc_comm_idx;
            cc_comm=spin_system.control.cc_comm;

            % Get the Hilbert dimension
            dim=size(drifts{1},1);

            % Loop over control sequence
            parfor n=1:(nsteps+1)

                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % First step is special
                if n==1

                    % Left pair of drifts
                    drift_pair={drifts{mod(n-1,ndrifts)+1},...
                                drifts{mod(n,ndrifts)+1}};

                    % Loop over controls
                    for k=1:nctrls

                        % Build the auxiliary matrix
                        [DL_first,~]=aux_mat(drift_pair,controls,cc_comm_idx,...
                                             cc_comm,dt(n),waveform(:,1),...
                                             waveform(:,2),k);

                        % Exponentiate the auxiliary matrix
                        auxmat=propagator(ss_parfor,DL_first,dt(n));

                        % Extract propagator and derivative
                        P_aux=auxmat(1:dim,1:dim);
                        dP=auxmat(1:dim,(1:dim)+dim);

                        % Apply the derivative to the density matrix
                        rho_deriv=dP*fwd_traj{n}*P_aux'+...
                                  P_aux*fwd_traj{n}*dP';

                        % Compute the derivative of the objective
                        grad_col(k)=hdot(bwd_traj{n+1},rho_deriv);

                    end

                % Last step is special
                elseif n==(nsteps+1)

                    % Right pair of drifts
                    drift_pair={drifts{mod(n-2,ndrifts)+1},...
                                drifts{mod(n-1,ndrifts)+1}};

                    % Loop over controls
                    for k=1:nctrls

                        % Build the auxiliary matrix
                        [~,DR_last]=aux_mat(drift_pair,controls,cc_comm_idx,...
                                            cc_comm,dt(n-1),...
                                            waveform(:,(end-1)),...
                                            waveform(:,end),k);

                        % Exponentiate the auxiliary matrix
                        auxmat=propagator(ss_parfor,DR_last,dt(n-1));

                        % Extract propagator and derivative
                        P_aux=auxmat(1:dim,1:dim);
                        dP=auxmat(1:dim,(1:dim)+dim);

                        % Apply the derivative to the density matrix
                        rho_deriv=dP*fwd_traj{n-1}*P_aux'+...
                                  P_aux*fwd_traj{n-1}*dP';

                        % Compute the derivative of the objective
                        grad_col(k)=hdot(bwd_traj{end},rho_deriv);

                    end

                % Middle steps
                else

                    % Loop over controls
                    for k=1:nctrls

                        % Left pair of drifts
                        drift_pair={drifts{mod(n-1,ndrifts)+1},...
                                    drifts{mod(n,ndrifts)+1}};

                        % Build the auxiliary matrix
                        [Right_DL,~]=aux_mat(drift_pair,controls,cc_comm_idx,...
                                             cc_comm,dt(n),waveform(:,n),...
                                             waveform(:,n+1),k);

                        % Exponentiate the auxiliary matrix
                        auxmat=propagator(ss_parfor,Right_DL,dt(n));

                        % Extract propagator and derivative
                        P_aux=auxmat(1:dim,1:dim);
                        dP=auxmat(1:dim,(1:dim)+dim);

                        % Apply the derivative to the density matrix
                        rho_deriv=dP*fwd_traj{n}*P_aux'+...
                                  P_aux*fwd_traj{n}*dP';

                        % Product rule: [dP2]*[P1]*rho part
                        grad_col(k)=grad_col(k)+hdot(bwd_traj{n+1},rho_deriv);

                        % Right pair of drifts
                        drift_pair={drifts{mod(n-2,ndrifts)+1},...
                                    drifts{mod(n-1,ndrifts)+1}};

                        % Build the auxiliary matrix
                        [~,Left_DR]=aux_mat(drift_pair,controls,cc_comm_idx,...
                                            cc_comm,dt(n-1),...
                                            waveform(:,n-1),...
                                            waveform(:,n),k);

                        % Exponentiate the auxiliary matrix
                        auxmat=propagator(ss_parfor,Left_DR,dt(n-1));

                        % Extract propagator and derivative
                        P_aux=auxmat(1:dim,1:dim);
                        dP=auxmat(1:dim,(1:dim)+dim);

                        % Apply the derivative to the density matrix
                        rho_deriv=dP*fwd_traj{n-1}*P_aux'+...
                                  P_aux*fwd_traj{n-1}*dP';

                        % Product rule: [P2]*[dP1]*rho part
                        grad_col(k)=grad_col(k)+hdot(bwd_traj{n},rho_deriv);

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

    % Check Hessian method
    if ~ismember(spin_system.control.method,{'newton','goodwin'})
        error('Hessian calculation methods are ''newton'' and ''goodwin''')
    end

    % Preallocate derivative propagators
    dP=cell(nctrls,nsteps);
    d2P=cell(nctrls,nctrls,nsteps);

    % Compute derivative propagators
    parfor n=1:nsteps

        % Preallocate local derivative arrays
        dP_col=cell(nctrls,1);
        d2P_block=cell(nctrls,nctrls);

        % Loop over control pairs
        for k=1:nctrls

            % First derivative of the propagator
            auxmat=dirdiff(ss_parfor,H{n},controls{k},dt(n),2);
            dP_col{k}=auxmat{2};

            % Second derivative of the propagator
            for j=1:nctrls
                auxmat=dirdiff(ss_parfor,H{n},{controls{k},controls{j}},dt(n),3);
                d2P_block{k,j}=auxmat{3};
            end

        end

        % Store derivative propagators
        dP(:,n)=dP_col;
        d2P(:,:,n)=d2P_block;

    end

    % Preallocate Hessian matrix
    hess=zeros(nctrls,nsteps,nctrls,nsteps,'like',1i);

    % Loop over derivative source intervals
    for m=1:nsteps

        % Preallocate local Hessian column
        hess_col=zeros(nctrls,nsteps,nctrls,1,'like',1i);

        % Loop over source controls
        for j=1:nctrls

            % Apply the first derivative at source interval
            rho_deriv=dP{j,m}*fwd_traj{m}*P{m}'+...
                      P{m}*fwd_traj{m}*dP{j,m}';

            % Loop over observer controls at the same interval
            for k=1:nctrls

                % Apply the second derivative at source interval
                rho_d2=d2P{k,j,m}*fwd_traj{m}*P{m}'+...
                       dP{k,m}*fwd_traj{m}*dP{j,m}'+...
                       dP{j,m}*fwd_traj{m}*dP{k,m}'+...
                       P{m}*fwd_traj{m}*d2P{k,j,m}';

                % Calculate non-mixed derivatives
                hess_col(k,m,j)=hdot(bwd_traj{m+1},rho_d2);

            end

            % Get the forward keyhole operator
            keyhole_forw=spin_system.control.keyholes{m};

            % Project the derivative trajectory
            if ~isempty(keyhole_forw), rho_deriv=keyhole_forw(rho_deriv); end

            % Loop over later observer intervals
            for n=(m+1):nsteps

                % Loop over observer controls
                for k=1:nctrls

                    % Apply the observer derivative
                    rho_mixed=dP{k,n}*rho_deriv*P{n}'+...
                              P{n}*rho_deriv*dP{k,n}';

                    % Calculate mixed derivatives
                    hess_col(k,n,j)=hdot(bwd_traj{n+1},rho_mixed);

                end

                % Propagate the derivative trajectory if needed
                if n<nsteps
                    rho_deriv=P{n}*rho_deriv*P{n}';
                    keyhole_forw=spin_system.control.keyholes{n};
                    if ~isempty(keyhole_forw), rho_deriv=keyhole_forw(rho_deriv); end
                end

            end

        end

        % Add to Hessian array
        hess(:,:,:,m)=hess_col;

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

% Return the trajectory (a huge array) only if needed
if (isfield(spin_system.control,'return_traj')&&spin_system.control.return_traj)||...
   any(ismember({'correlation_order','coherence_order',...
                 'local_each_spin',  'total_each_spin',...
                 'level_populations','trajectory'},spin_system.control.plotting(:)))
    traj_data.forward=fwd_traj;
else
    traj_data.forward=[];
end

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
function grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)
if ~ismember(spin_system.bas.formalism,{'zeeman-hilb'})
    error('this function requires a density matrix based formalism.');
end
if (~isnumeric(rho_init))||(~ismatrix(rho_init))||(size(rho_init,1)~=size(rho_init,2))
    error('rho_init must be a square matrix.');
end
if (~isnumeric(rho_targ))||(~ismatrix(rho_targ))||(size(rho_targ,1)~=size(rho_targ,2))
    error('rho_targ must be a square matrix.');
end
if (~ischar(fidelity_type))||(~ismember(fidelity_type,{'real','imag','square'}))
    error('fidelity_type must be ''real'', ''imag'', or ''square''.');
end
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
    if (size(drifts{n},1)~=size(rho_init,1))||...
       (size(drifts{n},1)~=size(rho_targ,1))
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

% Series of headlines in Le Moniteur Universel ("le journal de la pensée 
% officielle"), as Napoleon was approaching Paris between 9 and 21 March
% 1815, according to Alexandre Dumas:
%
%  "L'anthropophage est sorti de son repaire."
%
%  "L'ogre de Corse vient de débarquer au golfe Juan."
%
%  "Le tigre est arrivé à Gap."
%
%  "Le monstre a couché à Grenoble."
%
%  "Le tyran a traversé Lyon."
%
%  "L'usurpateur a été vu à soixante lieues de la capitale."
%
%  "Bonaparte s'avance à grands pas, mais il n'entrera jamais dans Paris."
%
%  "Napoléon sera demain sous nos remparts."
%
%  "L'empereur est arrivé à Fontainebleau."
%
%  "Sa Majesté Impériale et Royale a fait hier son entrée en son château 
%   des Tuileries au milieu de ses fidèles sujets."
%
% Dumas concludes: "This is the ultimate monument to journalism. It need not
% do anything else, for it won't do anything better."
