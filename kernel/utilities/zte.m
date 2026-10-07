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
% Note: sys.tols.zte_warr (default 1e-6) is an independent absolute
%       state error tolerance; it does not affect track selection.
%       The reported warranty bounds the 2-norm (and hence every
%       component) of exp(-1i*L*t)*rho-P*exp(-1i*P'*L*P*t)*P'*rho
%       for this supplied vector and fixed generator only. With
%       A=-1i*L, retained indices S, and discarded indices D, the
%       bound is exp(alpha*t)*(d+b*r*t), where d=norm(rho(D)),
%       r=norm(rho(S)), b>=norm(A(D,S),2), and alpha bounds the
%       non-negative logarithmic norm of A by Gershgorin discs.
%       Sparse norm scans and scalar bisection avoid eigensolvers.
%       This exact-arithmetic truncation bound excludes roundoff,
%       propagation error, future generators, other input vectors,
%       and frequency-domain spectra. A short warranty does not
%       imply that the actual error becomes large afterwards.
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
                            (nstates<1)||(mod(nstates,1)~=0)||(nstates>numel(rho)))
    error('nstates must be a positive integer not exceeding the state space dimension.');
end

% Initialise the unchanged-space warranty
warr_time=Inf;

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
     
    % Bound the initial discarded mass and retained-to-discarded leakage
    discarded=norm(rho(zero_track_mask),2);
    retained=norm(rho(~zero_track_mask),2);
    leakage=L(zero_track_mask,~zero_track_mask);
    leak_norm=sqrt(norm(leakage,1))*sqrt(norm(leakage,inf));

    % Distinguish an invalid initial projection from an exact invariant space
    if discarded>spin_system.tols.zte_warr
        warr_time=NaN;
    elseif (discarded>0)||((leak_norm>0)&&(retained>0))

        % Bound both full and compressed semigroups using the Hermitian part
        herm_part=(-1i*L)/2+(-1i*L)'/2;
        diagonal=real(diag(herm_part));
        herm_part=herm_part-spdiags(diagonal,0,size(L,1),size(L,2));
        growth=max(0,full(max(diagonal+sum(abs(herm_part),2))));

        % Solve the monotone bound in the log domain without exponent overflow
        log_leak=log(leak_norm)+log(retained);
        if (leak_norm==0)||(retained==0), log_leak=-Inf; end
        lower=0; upper=realmax;
        if (discarded>0)&&(growth>0)
            upper=min(upper,(log(spin_system.tols.zte_warr)-log(discarded))/growth);
        end
        if isfinite(log_leak)
            upper=min(upper,exp(log(spin_system.tols.zte_warr)-log_leak));
        end
        if isfinite(growth)&&(~isnan(log_leak))

            % Cover the floating-point exponent range and significand in bisection
            for k=1:(2-log2(realmin)+log2(realmax)-log2(eps))
                trial=lower+(upper-lower)/2;
                if (trial==lower)||(trial==upper), break; end
                log_terms=[log(discarded) log_leak+log(trial)];
                largest=max(log_terms);
                log_bound=growth*trial+largest+log(sum(exp(log_terms-largest)));
                if log_bound<=log(spin_system.tols.zte_warr)
                    lower=trial;
                else
                    upper=trial;
                end
            end
        end
        warr_time=lower;
        if (growth==0)&&(log_leak==-Inf), warr_time=Inf; end
    end

    % Report back to the user
    report(spin_system,['state space dimension reduced from ' num2str(size(L,2)) ...
                        ' to ' num2str(size(projector,2))]);
    
end

% Report the supplied-vector truncation warranty without certifying roundoff
if isnan(warr_time)
    report(spin_system,['ZTE warranty: no interval (not even 0 seconds) at tolerance ' ...
                        num2str(spin_system.tols.zte_warr,17) ...
                        '; initial discarded 2-norm exceeds tolerance.']);
else
    report(spin_system,['ZTE warranty: supplied-vector absolute 2-norm error <= ' ...
                        num2str(spin_system.tols.zte_warr,17) ' for 0 <= t <= ' ...
                        num2str(warr_time,17) ' seconds (fixed generator; excludes roundoff).']);
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


