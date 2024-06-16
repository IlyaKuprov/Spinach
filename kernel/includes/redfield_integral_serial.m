% Redfield integral evaluation, the serial path. This include is 
% called from within relaxation.m and follows the notation used
% in IK's paper:
%
%           http://dx.doi.org/10.1016/j.jmr.2010.12.004
%
% with the difference that the numerical quadrature method propo-
% sed there has been superceded by the much faster auxiliary mat-
% rix method described in:
%
%                http://dx.doi.org/10.1063/1.4928978
%
% This include is called when relaxation.m is not at the top of 
% the parallelisation call stack. When it is, the asynchronous
% include is called instead.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=redfield_integral_serial.m>

% Sum over spherical ranks
for n=1:numel(Q)
    
    % Sum over first projection index
    for k=1:(2*n+1)
        
        % Sum over second projection index
        for m=1:(2*n+1)
            
            % Only proceed if necessary
            if cheap_norm(Q{n}{k,m})>0
                
                % Sum over third projection index
                for p=1:(2*n+1)
                    
                    % Sum over fourth projection index
                    for q=1:(2*n+1)
                        
                        % Only proceed if necessary
                        if cheap_norm(Q{n}{p,q})>0
                            
                            % Get decay weights and rates for correlation function
                            [weights,rates,states]=corrfun(spin_system,n,k,m,p,q);
                            
                            % Loop over chemical species
                            for s=1:numel(states)
                                
                                % Loop over correlation function exponentials
                                for j=1:numel(weights{s})
                                    
                                    % Only proceed if necessary
                                    if abs(weights{s}(j))>0
                                        
                                        % Set the upper integration limit according to the accuracy goal
                                        upper_limit=-1.5*(1/rates{s}(j))*log(1/spin_system.tols.rlx_integration);
                    
                                        % Inform the user
                                        report(spin_system,['integrating L=' num2str(n) ', k=' num2str(n+1-k,'%+d') ... 
                                                            ', m=' num2str(n+1-m,'%+d') ', p=' num2str(n+1-p,'%+d') ...
                                                            ', q=' num2str(n+1-q,'%+d') ', chemical species ' num2str(s) '...']);
                                            
                                        % Kill the terms in L0 that are irrelevant on the time scale of the integration
                                        B=clean_up(spin_system,L0,spin_system.tols.rlx_integration/abs(upper_limit));
                                                            
                                        % Prepare the relevant matrices
                                        A=Q{n}{k,m}; C=Q{n}{p,q}'; D=B-1i*rates{s}(j)*speye(size(B));
                                        
                                        % Obliterate irrelevant elements
                                        A(~states{s},~states{s})=0; B(~states{s},~states{s})=0;
                                        C(~states{s},~states{s})=0; D(~states{s},~states{s})=0;
                                            
                                        % Compute the integral and clean up to the tolerance
                                        R_int=-weights{s}(j)*A*expmint(spin_system,B,C,D,upper_limit);
                                        R=clean_up(spin_system,R,1e-2*spin_system.tols.rlx_zero);
                                        
                                        % Add to the total
                                        R=R+R_int;
                                            
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
    
% Deallocate variables
clear('Q','L0','R_int','weights','rates','states');

% Report to the user
report(spin_system,'Redfield superoperator done.');

% Think not of what was asked, but of why. When you figure out why,
% you would know how to answer.
%
% Maxim Gorky

