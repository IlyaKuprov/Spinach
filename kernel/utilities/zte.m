% Zero track elimination function. Inspects the first few steps in the
% system trajectory and drops the states that did not get populated to 
% a user-specified tolerance. Syntax:
%
%                projector=zte(spin_system,L,rho,nstates)
%
% Parameters:
%
%      L       - the Liouvillian to be used for time 
%                propagation
%
%      rho     - the initial state to be used for 
%                time propagation
%
%      nstates - if this parameter is specified, only
%                nstates most populated states are kept,
%                irrespective of the tolerance parameter
%
% Output:
%
%      projector - projector matrix into the reduced space,
%                  to be used as follows: 
%
%                            L_reduced=P'*L*P
%                            rho_reduced=P'*rho;
%
% Note: default tolerance may be altered by setting sys.tols.zte_tol
%       variable before calling create.m 
%
% Note: further information on how this function works is available 
%       in IK's JMR paper on the subject
%
%               http://dx.doi.org/10.1016/j.jmr.2008.08.008
%
% Note: if tiny interactions or nearly equivalent spins are present,
%       it is best to disable zero track elimination by adding 'zte'
%       to the sys.disable cell array. 
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=zte.m>

function projector=zte(spin_system,L,rho,nstates)

% Validate the input
grumble(spin_system,L,rho);

% Run Zero Track Elimination
if ismember('zte',spin_system.sys.disable)
    
    % Skip if instructed to do so by the user
    report(spin_system,'WARNING - zero track elimination disabled, basis left unchanged.');
    
    % Return a unit matrix
    projector=1;

elseif nnz(rho)/numel(rho)>spin_system.tols.zte_maxden
    
    % Skip if the benefit is likely to be minor
    report(spin_system,'WARNING - too few zeros in the state vector, basis left unchanged.');
    
    % Return a unit matrix
    projector=1;
    
elseif norm(rho,1)<spin_system.tols.zte_tol
    
    % Skip if the state vector norm is too small for Krylov procedure
    report(spin_system,'WARNING - state vector norm below drop tolerance, basis left unchanged.');
    
    % Return a unit matrix
    projector=1;
    
else
    
    % Get the time step
    timestep=1/cheap_norm(L);
    
    % Do not allow infinite time step
    if isinf(timestep)
        report(spin_system,'zero Liouvillian supplied, using unit time step.'); timestep=1;
    end
    
    % Report to the user
    if exist('nstates','var')
        report(spin_system,['keeping ' num2str(nstates) ' states with the greatest trajectory weight.']); 
    else
        report(spin_system,['dropping states with amplitudes below ' num2str(spin_system.tols.zte_tol)...
                            ' within the first ' num2str(timestep*spin_system.tols.zte_nsteps) ' seconds.']);
    end
    report(spin_system,['a maximum of ' num2str(spin_system.tols.zte_nsteps) ...
                        ' steps shall be taken, ' num2str(timestep) ' seconds each.']);
    
    % Preallocate the trajectory
    trajectory=zeros(numel(rho),spin_system.tols.zte_nsteps,'like',1i);
    
    % Set the starting point
    trajectory(:,1)=rho;
    report(spin_system,['evolution step 0, active space dimension ' num2str(nnz(abs(trajectory(:,1))>spin_system.tols.zte_tol))]);
    
    % Compute trajectory steps with Krylov technique
    for n=2:spin_system.tols.zte_nsteps
        
        % Take a step forward
        trajectory(:,n)=step(spin_system,L,trajectory(:,n-1),timestep);
        
        % Analyze the trajectory
        prev_space_dim=nnz(max(abs(trajectory(:,1:(n-1))),[],2)>spin_system.tols.zte_tol);
        curr_space_dim=nnz(max(abs(trajectory),[],2)>spin_system.tols.zte_tol);
        
        % Inform the user
        report(spin_system,['evolution step ' num2str(n-1) ...
                            ', active space dimension ' num2str(curr_space_dim)]);
        
        % Terminate if done early
        if curr_space_dim==prev_space_dim, break; end
        
    end
    
    % Determine which tracks to drop
    if exist('nstates','var')
        
        % Determine state amplitudes
        amplitudes=max(abs(trajectory),[],2);
        
        % Sort the maximum amplitudes in descending order
        [~,index]=sort(amplitudes,'descend');
        
        % Drop all states beyond a given number
        zero_track_mask=true(size(rho));
        zero_track_mask(index(1:nstates))=false();
        
    else
        
        % Drop all states with maximum amplitude below the threshold 
        zero_track_mask=(max(abs(trajectory),[],2)<spin_system.tols.zte_tol);
        
    end
    
    % Take a unit matrix and delete the columns corresponding to zero tracks
    projector=speye(size(L)); projector(:,zero_track_mask)=[];
     
    % Report back to the user
    report(spin_system,['state space dimension reduced from ' num2str(size(L,2)) ...
                        ' to ' num2str(size(projector,2))]);
    
end

end

% Input validation function
function grumble(spin_system,L,rho)
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('zero track elimination is only available for zeeman-liouv and sphten-liouv formalisms.');
end
if (~isnumeric(L))||(~isnumeric(rho))
    error('both inputs must be numeric.');
end
if ~isvector(rho)
    error('single state vector expected, not a stack.');
end
if size(L,1)~=size(L,2)
    error('Liouvillian must be square.');
end
if size(L,2)~=size(rho,1)
    error('Liouvillian and state vector dimensions must be consistent.');
end
end

% Every great scientific truth goes through three stages. First, people say
% it conflicts with the Bible. Next they say it had been discovered before.
% Lastly they say they always believed it. 
%
% Louis Agassiz

