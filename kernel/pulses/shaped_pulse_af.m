% Shaped pulse in amplitude-frequency coordinates using Fokker-Planck
% formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jmr.2016.07.005).
% Syntax:
%
%  [rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
%                               rf_amp_list,rf_dur_list,rf_phi,...
%                               max_rank,method)
%
% Parameters:
%
%        L0          - drift Liouvillian that continues
%                      running in the background
%
%        Lx          - X projection of the RF operator
%
%        Ly          - Y projection of the RF operator
%
%        rho         - initial state vector or a horizontal 
%                      stack thereof
%
%        rf_frq_list - a vector of RF frequencies at each
%                      time slice (relative to the offsets
%                      and/or rotating frames that were
%                      used to make the background L0 that 
%                      you have supplied, Hz
%
%        rf_amp_list - a vector of RF amplitudes at each
%                      time slice, rad/s
%
%        rf_dur_list - a vector of time slice durations,
%                      in seconds
%
%        rf_phi      - RF phase of the first pulse slice
%
%        max_rank    - maximum rank of the Fokker-Planck
%                      theory, increase until the answer 
%                      stops changing, 2 is a good start
%
%        method      - propagation method, 'expv' for Krylov
%                      propagation, 'expm' for exponential
%                      propagation, 'evolution' for Spinach
%                      evolution function
%
% Outputs:
%
%        rho         - final state vector or a stack thereof
%
%        traj        - system trajectory as a [1 x (nsteps+1)] 
%                      cell array, the first point is the ini-
%                      tial condition
%
%        P           - effective pulse propagator (expensive),
%                      only available for the 'expm' method
%
% Note: the pulse is assumed to be piecewise-constant and should be 
%       supplied with sufficiently fine time discretisation to pro-
%       perly reproduce the waveform.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=shaped_pulse_af.m>

function [rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
                                      rf_amp_list,rf_dur_list,rf_phi,...
                                      max_rank,method)

% Set the defaults
if ~exist('method','var'), method='expv'; end

% Check consistency
grumble(L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,...
        rf_dur_list,rf_phi,max_rank,method);

% Get problem dimensions
spc_dim=2*max_rank+1; spn_dim=size(rho,1); stk_dim=size(rho,2);
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['state vector stack size         ' num2str(stk_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(spc_dim*spn_dim)]);

% Compute RF phases and Fourier derivative operator
[phases,d_dphi]=fourdif(spc_dim,1);

% Add the overall phase
phases=phases+rf_phi;

% Build the background Liouvillian
F0=polyadic({{opium(spc_dim,1),L0}});

% Build RF operators at all phases
F1=polyadic({{spdiags(cos(phases),0,spc_dim,spc_dim),Lx},...
             {spdiags(sin(phases),0,spc_dim,spc_dim),Ly}});

% Build RF phase turning generator
M=polyadic({{d_dphi,opium(spn_dim,1)}});

% Inflate polyadic representations
if ~ismember('polyadic',spin_system.sys.enable)
    F0=complex(inflate(F0)); 
    F1=complex(inflate(F1)); 
     M=complex(inflate(M));
end

% Project the state into the first block
proj=zeros(spc_dim,1); proj(1)=1; rho=kron(proj,rho);

% Upload pertinent things to GPU
if ismember('gpu',spin_system.sys.enable)
    F0=gpuArray(F0); F1=gpuArray(F1); 
    M=gpuArray(M);  rho=gpuArray(rho);
    location='GPU';
else
    location='CPU';
end

% Get the trajectory going
if nargout>1

    % Preallocate trajectory array
    traj=cell(1,numel(rf_dur_list)+1);

    % Fold the first point state stack back into Liouville space
    traj{1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2));

end

% Run the pulse
switch method
    
    case 'expv'
        
        % Start feedback timer
        feedback=tic(); 
        
        % Run Krylov propagation
        for n=1:numel(rf_frq_list)
            
            % Take a step forward
            rho=step(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,...
                                 rho,rf_dur_list(n));
            
            % Inform the user
            if (n==numel(rf_frq_list))||(toc(feedback)>1)
                report(spin_system,[location ' Krylov propagation step ' num2str(n) ...
                                    ' of ' num2str(numel(rf_frq_list)) ' done.']);
                feedback=tic();
            end

            % Store trajectory
            if nargout>1

                % Fold the state stack back into Liouville space
                traj{n+1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2)); 

            end 
            
        end

        % Refuse to compute the propagator
        if nargout>2, error('use ''expm'' method to get effective propagator.'); end
        
    case 'expm'
        
        % Get the total propagator going
        if nargout>2, P_tot=speye(size(F0)); end
                
        % Propagate the system
        for n=1:numel(rf_frq_list)

            % Compute the exponential propagator explicitly
            P=propagator(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,rf_dur_list(n));

            % Take the step
            rho=P*rho;

            % Store trajectory
            if nargout>1

                % Fold the state stack back into Liouville space
                traj{n+1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2));

            end 

            % Update pulse propagator in Fokker-Planck space
            if nargout>2
                P_tot=clean_up(spin_system,P*P_tot,spin_system.tols.prop_chop);
            end

        end
        
        % Project pulse propagator into Liouville space
        if nargout>2

            % Break the propagator up into spatial location blocks
            P_tot=mat2cell(P_tot,spn_dim*ones(1,spc_dim),...
                                 spn_dim*ones(1,spc_dim));

            % The effective propagator starts at first location and
            % collects the dynamics from all locations in the end
            P=P_tot{1,1};
            for n=2:spc_dim
                P=P+P_tot{n,1};
            end
            
        end
        
    case 'evolution'
        
        % Use the evolution function
        for n=1:numel(rf_frq_list)

            % Call the evolution function
            rho=evolution(spin_system,F0+rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,...
                                      [],rho,rf_dur_list(n),1,'final');

            % Store trajectory
            if nargout>1

                % Fold the state stack back into Liouville space
                traj{n+1}=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2)); 

            end

        end

        % Refuse to compute the propagator
        if nargout>2, error('use ''expm'' method to get effective propagator.'); end
        
    otherwise
        
        % Complain and bomb out
        error('unknown propagation method.');
        
end

% Fold the state stack back into Liouville space
rho=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2));

% Gather results from the GPU
if ismember('gpu',spin_system.sys.enable)
    rho=gather(rho);
    if nargout>1, traj=gather(traj); end
    if nargout>2, P=gather(P); end
end

end

% Consistency enforcement
function grumble(L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,...
                 rf_dur_list,rf_phi,max_rank,method)
if (~isnumeric(L0))||(~isnumeric(Lx))||(~isnumeric(Ly))||...
    (size(L0,1)~=size(L0,2))||(size(Lx,1)~=size(Lx,2))||...
    (size(Ly,1)~=size(Ly,2))||(size(L0,1)~=size(Lx,1))||...
    (size(Lx,1)~=size(Ly,1))
    error('L0, Lx, Ly must be square matrices of the same size.');
end
if size(L0,2)~=size(rho,1)
    error('matrix dimensions of L0 and rho must agree.');
end
if (~isnumeric(rf_frq_list))||(~isreal(rf_frq_list))||any(~isfinite(rf_frq_list))
    error('rf_frq_list must be a vector of real numbers.');
end
if (~isnumeric(rf_amp_list))||(~isreal(rf_amp_list))||any(~isfinite(rf_amp_list))
    error('rf_amp_list must be a vector of real numbers.');
end
if (~isnumeric(rf_dur_list))||(~isreal(rf_dur_list))||any(~isfinite(rf_dur_list))
    error('rf_amp_list must be a vector of real numbers.');
end
if (numel(rf_frq_list)~=numel(rf_amp_list))||...
   (numel(rf_amp_list)~=numel(rf_dur_list))
    error('rf_frq_list, rf_amp_list and rf_dur_list must have the same number of elements.');
end
if (~isnumeric(rf_phi))||(~isreal(rf_phi))||(~isfinite(rf_phi))||(numel(rf_phi)~=1)
    error('rf_phi must be a real number.');
end
if (~isnumeric(max_rank))||(~isreal(max_rank))||(~isfinite(max_rank))||...
   (numel(max_rank)~=1)||(max_rank<1)||mod(max_rank,1)
    error('max_rank must be a positive real integer.');
end
if ~ischar(method)
    error('method must be a character string.');
end
end

% A man gazing at the stars is at the mercy 
% of the puddles in the road.
%
% Alexander Smith

