% Time evolution function. Performs all types of time propagation with
% automatic trajectory level state space restriction. Syntax:
%
%         answer=evolution(spin_system,L,coil,rho,timestep,...
%                          nsteps,output,destination)
%
% Arguments for Liouville space calculations:
%
%      L      - the Liouvillian to be used during evolution
%
%      rho    - the initial state vector or a horizontal stack thereof
%
%      output - a string giving the type of evolution that is required
%
%                'final' - returns the final state vector or a horizontal
%                          stack thereof.
%
%                'trajectory' - returns the stack of state vectors giving
%                               the trajectory of the system starting from
%                               rho with the user-specified number of steps
%                               and step length.
%
%                'total'   - returns the integral of the observable trace
%                            from the simulation start to infinity. This
%                            option requires the presence of relaxation.
%
%                'refocus' - evolves the first vector for zero steps,
%                            second vector for one step, third vector for
%                            two steps, etc., consistent with the second
%                            stage of evolution in the indirect dimension
%                            after a refocusing pulse.
%
%                'observable' - returns the time dynamics of an observable
%                               as a vector (if starting from a single ini-
%                               tial state) or a matrix (if starting from a
%                               stack of initial states).
%
%                'multichannel' - returns the time dynamics of several
%                                 observables as rows of a matrix. Note
%                                 that destination state screening may be
%                                 less efficient when there are multiple
%                                 destinations to screen against.
%
%      coil   - the detection state, used when 'observable' is specified as
%               the output option. If 'multichannel' is selected, the coil
%               should contain multiple columns corresponding to individual
%               observable vectors.
%
%      destination - (optional) the state to be used for destination state
%                    screening.
%
% Arguments for Hilbert space calculations:
%
%       L         - Hamiltonian matrix
%
%       coil      - observable operator (if any)
%
%       rho       - initial density matrix
%
%       timestep  - duration of a single time step (seconds)
%
%       nsteps    - number of steps to take
%
%       output    - a string giving the type of evolution that is required
%
%                'final' - returns the final density matrix.
%
%                'trajectory' - returns a cell array of density matrices
%                               giving the trajectory of the system star-
%                               ting from rho with the user-specified num-
%                               ber of steps and step length.
%
%                'refocus' - evolves the first matrix for zero steps,
%                            second matrix for one step, third matrix for
%                            two steps, etc., consistent with the second
%                            stage of evolution in the indirect dimension
%                            after a refocusing pulse.
%
%                'observable' - returns the time dynamics of an observable
%                               as a vector.
%
%       destination - this argument is ignored.
%
% Outputs:
% 
%       answer - a vector, a matrix, or a cell array or matrices,
%                depending on the options set during the call
%
% Calculation of final states and observables in Hilbert space is parallel-
% ized and tested all the way to 128-core (16 nodes, 8 cores each) configu-
% rations. Parallelization of the trajectory calculation does not appear to
% yield any benefits due to large amount of inter-thread communication. See
% http://dx.doi.org/10.1063/1.3679656 for further information.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% ohad.levinkron@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=evolution.m>

function answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)

% Check consistency
grumble(L,coil,rho,timestep,nsteps,output);

% Call Krylov propagation for polyadics
if isa(L,'polyadic')
    report(spin_system,'polyadic generator received, forwarding to krylov()...');
    answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output); return;
end

% Gather state vectors from GPUs
if isa(rho,'gpuArray'), rho=gather(rho); end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case {'sphten-liouv','zeeman-liouv'}
        
        % Apply trajectory-level reduction algorithms
        report(spin_system,'trying to reduce the problem dimension...');
        
        % Decide the screening algorithm
        if ismember('dss',spin_system.sys.disable)
            
            % If DSS is disabled, run forward reduction
            report(spin_system,'WARNING - destination state screening is disabled.');
            projectors=reduce(spin_system,L,rho);
            
        elseif exist('destination','var')&&(~isempty(destination))
            
            % If destination state is supplied, use it for DSS
            report(spin_system,'destination state screening using supplied state.');
            projectors=reduce(spin_system,L',destination);
            
        elseif ismember(output,{'observable'})
            
            % If coil state is supplied, use it for DSS
            report(spin_system,'destination state screening using coil state.');
            projectors=reduce(spin_system,L',coil);
            
        else
            
            % Default to the usual forward screening
            report(spin_system,'destination state screening is not applicable.');
            projectors=reduce(spin_system,L,rho);
            
        end
        
        % Count subspaces
        nsubs=numel(projectors);
        
        % Parallel strategy
        if nsubs>1
            nworkers=poolsize;
        else
            nworkers=0;
        end
        
        % Run the evolution
        switch output
            
            case 'final'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); rho_sub=cell(nsubs,1); 
                rows=cell(nsubs,1); cols=cell(nsubs,1); vals=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop in parallel over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    
                    % Decide if Krylov should be used
                    if ((size(L_loc,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Ignore user timing and use optimal one
                        nsteps_opt=ceil(log2(size(L_loc,1)/(size(rho_loc,2)*log(2)))); 
                        nsteps_opt=max([1 nsteps_opt]); total_time=timestep*nsteps;
                        timestep_loc=total_time/nsteps_opt;
                        
                        % Inform the user
                        report(spin_system,['dim(L)=' num2str(size(L_loc,1)) ...
                                            ', rho stack size ' num2str(size(rho_loc,2)) ', ' ...
                                            num2str(nsteps_opt) ' substeps in optimal stepping.']); 
                  
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep_loc);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move things to GPU
                            P=gpuArray(P); rho_loc=gpuArray(rho_loc);
                            
                            % Run the propagation
                            for n=1:nsteps_opt
                                rho_loc=P*rho_loc;
                            end
                            
                            % Gather the result
                            rho_loc=gather(rho_loc);
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Run the propagation
                            for n=1:nsteps_opt
                                rho_loc=P*rho_loc;
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_loc=krylov(spin_system,L_loc,[],rho_loc,timestep,nsteps,'final');
                        
                    end
                    
                    % Project back into the full space and convert into i,j,v
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    [rows{sub},cols{sub},vals{sub}]=find(proj_loc*sparse(rho_loc));
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; %#ok<NASGU>
                    
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
                % Deallocate memory
                clear('L_sub','rho_sub','projectors');
                
                % Recombine the subspaces
                report(spin_system,'recombining subspaces...');
                rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals);
                answer=sparse(rows,cols,vals,size(rho,1),size(rho,2));
                answer=clean_up(spin_system,answer,spin_system.tols.zte_tol);
                report(spin_system,'data retrieval finished.');
                
            case 'total'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=full(projectors{sub}'*rho);
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Start with zero
                answer=0;
                
                % Loop in parallel over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; coil_loc=coil_sub{sub}; 
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    coil_loc=clean_up(spin_system,coil_loc,spin_system.tols.zte_tol);
                    
                    % Compute the integral
                    report(spin_system,'integrating the trajectory...');
                    answer_loc=real(coil_loc'*((1i*L_loc)\rho_loc));
                    
                    % Add the subspace to the total
                    answer=answer+answer_loc;
                    
                    % Update the user
                    report(spin_system,'integration finished.');
                    
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            case 'trajectory'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); rho_sub=cell(nsubs,1);
                rows=cell(nsubs,1); cols=cell(nsubs,1); vals=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop in parallel over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                     % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    
                    % Propagate the system
                    if ((size(L_loc,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Preallocate the answer
                        answer=zeros([size(rho_loc,1) (nsteps+1)],'like',1i);
                  
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                    
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move things to GPU
                            rho_loc=gpuArray(rho_loc); P=gpuArray(P);
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer(:,n)=gather(rho_loc);
                                rho_loc=P*rho_loc;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer(:,n)=rho_loc;
                                rho_loc=P*rho_loc;
                            end
                            
                        end

                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer=krylov(spin_system,L_loc,[],rho_loc,timestep,nsteps,'trajectory')
                        
                    end
                    
                    % Project back into the full space
                    answer=clean_up(spin_system,answer,spin_system.tols.zte_tol);
                    [rows{sub},cols{sub},vals{sub}]=find(proj_loc*sparse(answer));
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; answer=[]; %#ok<NASGU>
                    
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
                % Deallocate memory
                clear('L_sub','rho_sub','projectors');
                
                % Recombine the subspaces
                report(spin_system,'recombining subspaces...');
                rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals);
                answer=sparse(rows,cols,vals,size(rho,1),nsteps+1);
                answer=clean_up(spin_system,answer,spin_system.tols.zte_tol);
                report(spin_system,'data retrieval finished.');
                
            case 'refocus'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); rho_sub=cell(nsubs,1);
                rows=cell(nsubs,1); cols=cell(nsubs,1); vals=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    
                    % Propagate the system
                    if ((size(L_loc,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move propagator to GPU
                            P=gpuArray(P);
                            
                            % Loop over stack elements
                            for n=2:size(rho_loc,2)               % done this way because sparse gpuArrays
                                                                  % cannot be indexed into in R2017b - IK
                                % Get the current state vector
                                rho_cur=rho_loc(:,n);
                                
                                % Move it to GPU
                                rho_cur=gpuArray(rho_cur);
                                
                                % Propagate it
                                for k=2:n
                                    rho_cur=P*rho_cur;
                                end
                                
                                % Gather and reassign
                                rho_loc(:,n)=gather(rho_cur);
                                
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=2:size(rho_loc,2)
                                rho_loc(:,n:end)=P*rho_loc(:,n:end);
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_loc=krylov(spin_system,L_loc,[],rho_loc,timestep,[],'refocus');
                        
                    end
                    
                    % Project back into the full space
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    [rows{sub},cols{sub},vals{sub}]=find(proj_loc*sparse(rho_loc));
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; %#ok<NASGU>
                    
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
                % Deallocate memory
                clear('L_sub','rho_sub','projectors');
                
                % Recombine the subspaces
                report(spin_system,'recombining subspaces...');
                rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals);
                answer=sparse(rows,cols,vals,size(rho,1),nsteps+1);
                answer=clean_up(spin_system,answer,spin_system.tols.zte_tol);
                report(spin_system,'data retrieval finished.');
                
            case 'observable'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Preallocate the answer
                answer=zeros([(nsteps+1) size(rho,2)],'like',1i);
                
                % Loop over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; coil_loc=coil_sub{sub}; 
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    coil_loc=clean_up(spin_system,coil_loc,spin_system.tols.zte_tol);
                    
                    % Propagate the system
                    if ((size(L_loc,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Preallocate the local answer
                        answer_loc=zeros([(nsteps+1) size(rho_loc,2)],'like',1i);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Upload things to GPU
                            P=gpuArray(P); coil_loc=gpuArray(coil_loc); rho_loc=gpuArray(rho_loc);
                            
                            % Propagate and detect
                            for n=1:(nsteps+1)
                                answer_loc(n,:)=gather(coil_loc'*rho_loc);
                                rho_loc=P*rho_loc;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate and detect
                            for n=1:(nsteps+1)
                                answer_loc(n,:)=coil_loc'*rho_loc;
                                rho_loc=P*rho_loc;
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_loc=krylov(spin_system,L_loc,coil_loc,rho_loc,timestep,nsteps,'observable');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_loc;
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; P=[]; answer_loc=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
            case 'multichannel'
                
                % Create arrays of projections
                L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Preallocate the answer
                answer=zeros([size(coil,2) (nsteps+1)],'like',1i);
                
                % Loop over independent subspaces
                parfor (sub=1:nsubs,nworkers)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(numel(projectors)) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; coil_loc=coil_sub{sub}; 
                    L_loc=clean_up(spin_system,L_loc,spin_system.tols.liouv_zero);
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    coil_loc=clean_up(spin_system,coil_loc,spin_system.tols.zte_tol);
                    
                    % Propagate the system
                    if ((size(L_loc,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Preallocate the local answer
                        answer_loc=zeros([size(coil_loc,2) (nsteps+1)],'like',1i);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Upload things to GPU
                            P=gpuArray(P); rho_loc=gpuArray(rho_loc); coil_loc=gpuArray(coil_loc);
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_loc(:,n)=gather(coil_loc'*rho_loc);
                                rho_loc=P*rho_loc;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_loc(:,n)=coil_loc'*rho_loc;
                                rho_loc=P*rho_loc;
                            end
                        
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_loc=krylov(spin_system,L_loc,coil_loc,rho_loc,timestep,nsteps,'multichannel');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_loc;
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; P=[]; answer_loc=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
            otherwise
                
                % Complain and bomb out
                error('invalid output option.');
                
        end
        
    case 'zeeman-hilb'
        
        % Run the reduction
        projectors=reduce(spin_system,L,rho);
        
        % Run the evolution
        switch output
            
            case 'final'
                
                % Ignore user timing
                total_time=timestep*nsteps;
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,total_time);
                        
                        % Report to the user
                        report(spin_system,'propagating the system...');
                        
                        % Parallelise if possible
                        if iscell(rho)
                            
                            % Run in parallel over cells
                            parfor k=1:numel(rho)
                                rho{k}=P*rho{k}*P'
                            end
                            
                            % Return the result
                            answer=rho;
                            
                        else
                            
                            % Just apply the propagator
                            answer=P*rho*P';
                            
                        end
                        
                        % Report to the user
                        report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            case 'observable'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % If a stack is received, run 2D acquisition
                        if iscell(rho)
                            
                            % Compute the exponential propagator
                            P=propagator(spin_system,L,timestep);
                            
                            % Preallocate the answer
                            answer=zeros([(nsteps+1) numel(rho)],'like',1i);
                            
                            % Report to the user
                            report(spin_system,'propagating the system (parfor over stack)...');
                            
                            % Loop over the elements of the stack
                            parfor k=1:numel(rho)
                                
                                % Grab local copies
                                rho_local=rho{k};
                                ans_local=zeros([(nsteps+1) 1],'like',1i);
                                
                                % Loop over time steps
                                for n=1:(nsteps+1)
                                
                                    % Compute the observable
                                    ans_local(n)=hdot(coil,rho_local);
                                
                                    % Step forward
                                    rho_local=P*rho_local*P';
                                    
                                end
                                
                                % Assign the slice
                                answer(:,k)=ans_local;
                                
                            end
                        
                        % If inside a parallel loop, avoid spmd
                        elseif isworkernode
                            
                            % Compute the exponential propagator
                            P=propagator(spin_system,L,timestep);
                            
                            % Preallocate the answer
                            answer=zeros([(nsteps+1) 1],'like',1i);
                            
                            % Report to the user
                            report(spin_system,'propagating the system (worker node)...');
                            
                            % Loop over time steps
                            for n=1:(nsteps+1)
                                
                                % Compute the observable
                                answer(n)=hdot(coil,rho);
                                
                                % Step forward
                                rho=P*rho*P';
                                
                            end
                            
                        % Otherwise use Kuprov-Edwards parallel split
                        else
                            
                            % Preallocate the answer
                            answer=zeros([(nsteps+1) 1],'like',1i);
                            
                            % Report to the user
                            report(spin_system,'propagating the system (Kuprov-Edwards parallel split)...');
                            
                            % Loop over the irreps
                            for sub=1:numel(projectors)
                                
                                % Report to the user
                                report(spin_system,['irreducible representation ' int2str(sub) ...
                                                    '/' int2str(numel(projectors)) '...']);
                                
                                % Project the calculation
                                L_sub=projectors{sub}'*L*projectors{sub};
                                rho_sub=projectors{sub}'*rho*projectors{sub};
                                coil_sub=projectors{sub}'*coil*projectors{sub};
                            
                                % Compute the exponential propagator
                                P_sub=propagator(spin_system,L_sub,timestep);
                                
                                % Distribute the initial condition
                                rho_sub=distrib_dim(rho_sub,2);
                            
                                % Distribute covector array
                                cov_sub=speye(size(rho_sub));
                                cov_sub=distrib_dim(cov_sub,2);
                            
                                % Parallel processing
                                spmd
                                
                                    % Localize the problem and save memory
                                    rho_local=getLocalPart(rho_sub); rho_sub=[]; %#ok<NASGU>
                                    cov_local=getLocalPart(cov_sub); cov_sub=[]; %#ok<NASGU>
                                    fid_local=zeros([(nsteps+1) 1],'like',1i);
                                
                                    % Loop over time steps
                                    for n=1:(nsteps+1)
                                    
                                        % Write local fid
                                        fid_local(n)=hdot(coil_sub*cov_local,rho_local);
                                    
                                        % Step forward on rho and cov
                                        rho_local=P_sub*rho_local;
                                        cov_local=P_sub*cov_local;
                                    
                                    end
                                
                                    % Collect the results
                                    answer_sub=spmdPlus(fid_local,1);
                                
                                end
                            
                                % Return the result
                                answer=answer+answer_sub{1};
                                
                            end
                        
                        end
                        
                        % Report to the user
                        report(spin_system,'propagation finished.');
                                                     
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            case 'trajectory'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Preallocate the answer
                        answer=cell(1,nsteps+1);
                        
                        % Compute the trajectory
                        report(spin_system,'propagating the system...');
                        for n=1:(nsteps+1)
                            answer{n}=rho; rho=P*rho*P';
                        end
                        report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            case 'refocus'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Propagate the system
                        report(spin_system,'propagating the system...');
                        parfor k=2:numel(rho)
                            
                            % Grab a local copy
                            rho_local=rho{k};
                            
                            % Propagate local copy
                            for n=2:k
                                rho_local=P*rho_local*P';
                            end
                            
                            % Return local copy
                            rho{k}=rho_local;
                            
                        end
                        answer=rho; report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            otherwise
                
                % Complain and bomb out
                error('unknown evolution option for zeeman-hilb formalism.');
                
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% Consistency enforcement
function grumble(L,coil,rho,timestep,nsteps,output)
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

% Degrees of ability vary, but the basic principle remains the same: the
% degree of a man's independence, initiative and personal love for his work
% determines his talent as a worker and his worth as a man. Independence is
% the only gauge of human virtue and value. What a man is and makes of him-
% self; not what he has or hasn't done for others. There is no substitute
% for personal dignity. There is no standard of personal dignity except
% independence.
%
% Ayn Rand, "The Fountainhead"

