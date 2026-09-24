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
%      rho     - initial state column or horizontal stack of
%                state columns to be used for time propagation
%
%      nstates - if this parameter is specified, only
%                nstates most populated states are kept,
%                irrespective of the tolerance parameter; stacks
%                are ranked by maximum amplitude over columns
%                and sampled times
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

% Validate the number of states if it is specified
if exist('nstates','var')&&((~isnumeric(nstates))||(~isreal(nstates))||(~isscalar(nstates))||...
                            (nstates<1)||(mod(nstates,1)~=0)||(nstates>size(rho,1)))
    error('nstates must be a positive integer not exceeding the state space dimension.');
end

% Run Zero Track Elimination
if ismember('zte',spin_system.sys.disable)
    
    % Skip if instructed to do so by the user
    report(spin_system,'WARNING - zero track elimination disabled, basis left unchanged.');
    
    % Return a unit matrix
    projector=1;

elseif nnz(any(rho,2))/size(rho,1)>spin_system.tols.zte_maxden
    
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
    
    % Track the maximum amplitude of each actual state column
    amplitudes=abs(rho);
    report(spin_system,['evolution step 0, active space dimension ' ...
                        num2str(nnz(any(amplitudes>spin_system.tols.zte_tol,2)))]);
    
    % Compute trajectory steps with Krylov technique
    for n=2:spin_system.tols.zte_nsteps
        
        % Record each column's active dimension before propagation
        prev_space_dim=sum(amplitudes>spin_system.tols.zte_tol,1);

        % Take a step forward with the true state columns
        rho=step(spin_system,L,rho,timestep);

        % Analyse the trajectories without mixing phases or columns
        amplitudes=max(amplitudes,abs(rho));
        curr_space_dim=sum(amplitudes>spin_system.tols.zte_tol,1);
        
        % Inform the user
        report(spin_system,['evolution step ' num2str(n-1) ...
                            ', active space dimension ' ...
                            num2str(nnz(any(amplitudes>spin_system.tols.zte_tol,2)))]);
        
        % Terminate if done early
        if all(curr_space_dim==prev_space_dim), break; end
        
    end
    
    % Screen the union of the actual column trajectories
    amplitudes=max(amplitudes,[],2);

    % Determine which tracks to drop
    if exist('nstates','var')
        
        % Sort the maximum amplitudes in descending order
        [~,index]=sort(amplitudes,'descend');
        
        % Drop all states beyond a given number
        zero_track_mask=true(size(rho,1),1);
        zero_track_mask(index(1:nstates))=false();
        
    else
        
        % Drop all states with maximum amplitude below the threshold 
        zero_track_mask=(amplitudes<spin_system.tols.zte_tol);
        
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
if (~ismatrix(rho))||(size(rho,2)==0)
    error('state column or horizontal stack of state columns expected.');
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

