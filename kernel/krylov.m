% Krylov propagation function. Avoids matrix exponentiation, but can be 
% slow. Should be used when the Liouvillian exponential does not fit in-
% to the system memory, but the Liouvillian itself does. Syntax:
%
%     answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)
%
% Arguments:
%
%   L      - the Liouvillian to be used during evolution
%
%   rho    - the initial state vector or a horizontal stack thereof
%
%   output - a string giving the type of evolution that is required
%
%             'final' - returns the final state vector or a horizontal
%                       stack thereof.
%
%             'trajectory' - returns the stack of state vectors giving
%                            the trajectory of the system starting from
%                            rho with the user-specified number of steps
%                            and step length.
%
%             'total'   - returns the integral of the observable trace
%                         from the simulation start to infinity. This
%                         option requires the presence of relaxation.
%
%             'refocus' - evolves the first vector for zero steps,
%                         second vector for one step, third vector for
%                         two steps, etc., consistent with the second
%                         stage of evolution in the indirect dimension
%                         after a refocusing pulse.
%
%             'observable' - returns the time dynamics of an observable
%                            as a vector (if starting from a single ini-
%                            tial state) or a matrix (if starting from a
%                            stack of initial states).
%
%             'multichannel' - returns the time dynamics of several
%                              observables as rows of a matrix. Note
%                              that destination state screening may be
%                              less efficient when there are multiple
%                              destinations to screen against.
%
%   coil   - the detection state, used when 'observable' is specified as
%            the output option. If 'multichannel' is selected, the coil
%            should contain multiple columns corresponding to individual
%            observable vectors.
%
% Outputs:
%
%      answer - a vector or a matrix, depending on the options set during
%               the call.
%
% Note: this function does not support Hilbert space formalisms.
%
% Note: we initially had a faithful implementation of the Krylov process
%       here - subspace, orthogonalisation, projection, etc., but in all
%       our testing it was much inferior to the reordered Taylor process
%       that is currently implemented below.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=krylov.m>

function answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)

% Check consistency
grumble(spin_system,L,coil,rho,timestep,nsteps,output);

% Upload data to GPU or optimise layout
if ismember('gpu',spin_system.sys.enable)
    L=gpuArray(L); rho=gpuArray(full(rho));
    coil=gpuArray(coil); location='GPU';
else
    location='CPU'; rho=full(rho);
end

% Start feedback timer
feedback=tic();

% Decide the output type
switch output
    
    case 'final'
        
        % Take the step, ignoring user settings
        rho=step(spin_system,L,rho,timestep*nsteps);
            
        % Assign the answer
        answer=gather(rho);
                
    case 'trajectory'
        
        % Preallocate the answer and set the starting point
        answer=zeros([size(rho,1) (nsteps+1)],'like',1i);
        answer(:,1)=gather(rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Take next step
            rho=step(spin_system,L,rho,timestep);
            
            % Assign the answer
            answer(:,n+1)=gather(rho);
            
            % Inform the user
            if (n==nsteps)||(toc(feedback)>1)
                report(spin_system,[location ' Krylov step ' num2str(n)...
                                    ' out of ' num2str(nsteps) ' done.']);
                feedback=tic();
            end
            
        end
        
    case 'refocus'
        
        % Number of steps is fixed
        nsteps=size(rho,2);
        
        % Loop over steps
        for n=2:nsteps

            % Take next step
            rho(:,n:end)=step(spin_system,L,rho(:,n:end),timestep);
            
            % Inform the user
            if (n==nsteps)||(toc(feedback)>1)
                report(spin_system,[location ' Krylov step ' num2str(n)...
                                    ' out of ' num2str(nsteps) ' done.']);
                feedback=tic();
            end
            
        end
        
        % Assign the answer
        answer=gather(rho);
        
    case 'observable'
        
        % Preallocate the answer
        answer=zeros([(nsteps+1) size(rho,2)],'like',1i);
  
        % Set the initial point
        answer(1,:)=gather(coil'*rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Take next step
            rho=step(spin_system,L,rho,timestep);
            
            % Assign the answer
            answer(n+1,:)=gather(coil'*rho);
            
            % Inform the user
            if (n==nsteps)||(toc(feedback)>1)
                report(spin_system,[location ' Krylov step ' num2str(n)...
                                    ' out of ' num2str(nsteps) ' done.']);
                feedback=tic();
            end
            
        end
        
    case 'multichannel'
        
        % Preallocate the answer
        answer=zeros([size(coil,2) (nsteps+1)],'like',1i);
        
        % Set the initial point
        answer(:,1)=answer(:,1)+gather(coil'*rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Take next step
            rho=step(spin_system,L,rho,timestep);
            
            % Assign the answer
            answer(:,n+1)=gather(coil'*rho);
            
            % Inform the user
            if (n==nsteps)||(toc(feedback)>1)
                report(spin_system,[location ' Krylov step ' num2str(n)...
                                    ' out of ' num2str(nsteps) ' done.']);
                feedback=tic();
            end
            
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown output option.');
        
end

end

% Consistency enforcement
function grumble(spin_system,L,coil,rho,timestep,nsteps,output)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function only works in Lioville space.');
end
if ~isnumeric(L)
    error('Liouvillian must be numeric.');
end
if ~isnumeric(coil)
    error('coil argument must be numeric.');
end
if (~isnumeric(rho))&&(~iscell(rho))
    error('rho argument must either be numeric or a cell array');
end
if ~isnumeric(timestep)
    error('timestep argument must be numeric.');
end
if ~isnumeric(nsteps)
    error('nsteps must be numeric.');
end
if (~ischar(output))||(~ismember(output,{'observable','final',...
   'trajectory','total','multichannel','refocus'}))
    error('observable argument must be a valid character string.');
end
end

% Freedom can be achieved and retained only by sober men who
% take humanity as it is, not as humanity should be.
%
% Russell Kirk

