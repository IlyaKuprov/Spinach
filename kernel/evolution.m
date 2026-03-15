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
% ilya.kuprov@weizmann.ac.il
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

        % Run the evolution
        switch output
            
            case 'final'
                
                % Preallocate the answer
                answer=zeros(size(rho),'like',1i);
                                
                % Over subspaces
                for sub=1:nsubs

                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                    % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                    
                    % Decide if Krylov should be used instead
                    if ((size(L_sub,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Ignore user-specified timing and use the optimal one
                        nsteps_opt=ceil(log2(size(L_sub,1)/(size(rho_sub,2)*log(2)))); 
                        nsteps_opt=max([1 nsteps_opt]); total_time=timestep*nsteps;
                        timestep_opt=total_time/nsteps_opt;
                        
                        % Inform the user
                        report(spin_system,['dim(L)=' num2str(size(L_sub,1))             ...
                                            ', rho stack size ' num2str(size(rho_sub,2)) ...
                                            ', ' num2str(nsteps_opt)                     ...
                                            ' substeps in optimal stepping.']); 
                  
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_sub,timestep_opt);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move things to GPU
                            rho_sub=gpuArray(rho_sub);
                            P=gpuArray(P);
                            
                            % Run the propagation
                            for n=1:nsteps_opt
                                rho_sub=P*rho_sub;
                            end
                            
                            % Gather the result
                            rho_sub=gather(rho_sub);
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Run the propagation
                            for n=1:nsteps_opt
                                rho_sub=P*rho_sub;
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_sub=krylov(spin_system,L_sub,[],rho_sub,timestep,nsteps,'final');
                        
                    end
                    
                    % Project back and add to the answer
                    answer=answer+projectors{sub}*rho_sub;
                    
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
            case 'total'
                
                % Start
                answer=0;
                
                % Over subspaces
                for sub=1:nsubs

                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    coil_sub=full(projectors{sub}'*coil);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                    % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                   
                    % Compute the integral and add it to the answer
                    report(spin_system,'integrating the trajectory...');
                    answer=answer+1i*coil_sub'*(L_sub\rho_sub);
                    
                    % Update the user
                    report(spin_system,'integration finished.');
                    
                end
                
            case 'trajectory'
                
                % Preallocate the answer
                answer=zeros(size(rho,1),nsteps+1,'like',1i);

                % Over subspaces
                for sub=1:nsubs

                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                    % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                    
                    % Decide if Krylov should be used instead
                    if ((size(L_sub,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Preallocate the answer in the current subspace
                        answer_subs=zeros([size(rho_sub,1) (nsteps+1)],'like',1i);
                  
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_sub,timestep);
                    
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move things to GPU
                            rho_sub=gpuArray(rho_sub); P=gpuArray(P);
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_subs(:,n)=gather(rho_sub);
                                rho_sub=P*rho_sub;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_subs(:,n)=rho_sub;
                                rho_sub=P*rho_sub;
                            end
                            
                        end

                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_subs=krylov(spin_system,L_sub,[],rho_sub,timestep,nsteps,'trajectory');
                        
                    end
                    
                    % Project back and add to the answer
                    answer=answer+projectors{sub}*answer_subs;
                    
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
            case 'refocus'
                
                % Preallocate the answer
                answer=zeros(size(rho),'like',1i);
                
                % Over subspaces
                for sub=1:nsubs

                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                    % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                    
                    % Decide if Krylov should be used instead
                    if ((size(L_sub,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_sub,timestep);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Move propagator to GPU
                            P=gpuArray(P);
                            
                            % Loop over stack elements
                            for n=2:size(rho_sub,2)               % done this way because sparse gpuArrays
                                                                  % cannot be indexed into in R2017b - IK
                                % Get the current state vector
                                rho_cur=rho_sub(:,n);
                                
                                % Move it to GPU
                                rho_cur=gpuArray(rho_cur);
                                
                                % Propagate it
                                for k=2:n
                                    rho_cur=P*rho_cur;
                                end
                                
                                % Gather and reassign
                                rho_sub(:,n)=gather(rho_cur);
                                
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=2:size(rho_sub,2)
                                rho_sub(:,n:end)=P*rho_sub(:,n:end);
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_sub=krylov(spin_system,L_sub,[],rho_sub,timestep,[],'refocus');
                        
                    end
                    
                    % Project back and add to the answer
                    answer=answer+projectors{sub}*rho_sub;
                   
                end
                
                % Update the user
                report(spin_system,'propagation finished.');
                
            case 'observable'
                
                % Preallocate the answer
                answer=zeros([(nsteps+1) size(rho,2)],'like',1i);
                
                % Loop over independent subspaces
                for sub=1:nsubs

                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    coil_sub=full(projectors{sub}'*coil);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                     % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                    
                    % Propagate the system
                    if ((size(L_sub,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_sub,timestep);
                        
                        % Preallocate the local answer
                        answer_sub=zeros([(nsteps+1) size(rho_sub,2)],'like',1i);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Upload things to GPU
                            coil_sub=gpuArray(coil_sub); 
                            rho_sub=gpuArray(rho_sub);
                            P=gpuArray(P);
                            
                            % Propagate and detect
                            for n=1:(nsteps+1)
                                answer_sub(n,:)=gather(coil_sub'*rho_sub);
                                rho_sub=P*rho_sub;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate and detect
                            for n=1:(nsteps+1)
                                answer_sub(n,:)=coil_sub'*rho_sub;
                                rho_sub=P*rho_sub;
                            end
                            
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_sub=krylov(spin_system,L_sub,coil_sub,rho_sub,timestep,nsteps,'observable');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_sub;
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
            case 'multichannel'
                
                % Preallocate the answer
                answer=zeros([size(coil,2) (nsteps+1)],'like',1i);
                
                % Over subspaces
                for sub=1:nsubs
                    
                    % Project into the current subspace
                    L_sub=projectors{sub}'*L*projectors{sub};
                    rho_sub=full(projectors{sub}'*rho);
                    coil_sub=full(projectors{sub}'*coil);
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ...
                                        ' of ' num2str(nsubs) '...']);
                    
                     % Clean up the projection of the evolution generator
                    L_sub=clean_up(spin_system,L_sub,spin_system.tols.liouv_zero);
                    
                    % Propagate the system
                    if ((size(L_sub,2)<spin_system.tols.krylov_tol)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_sub,timestep);
                        
                        % Preallocate the local answer
                        answer_sub=zeros([size(coil_sub,2) (nsteps+1)],'like',1i);
                        
                        % Adapt to the target device
                        if ismember('gpu',spin_system.sys.enable)
                            
                            % Inform the user
                            report(spin_system,'propagating the system on GPU...');
                            
                            % Upload things to GPU
                            P=gpuArray(P); 
                            rho_sub=gpuArray(rho_sub); 
                            coil_sub=gpuArray(coil_sub);
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_sub(:,n)=gather(coil_sub'*rho_sub);
                                rho_sub=P*rho_sub;
                            end
                            
                        else
                            
                            % Inform the user
                            report(spin_system,'propagating the system on CPU...');
                            
                            % Propagate the system
                            for n=1:(nsteps+1)
                                answer_sub(:,n)=coil_sub'*rho_sub;
                                rho_sub=P*rho_sub;
                            end
                        
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_sub=krylov(spin_system,L_sub,coil_sub,rho_sub,timestep,nsteps,'multichannel');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_sub;
                    
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
                                
                                % Partition the density matrix columns into
                                % contiguous segments for parallel processing.
                                ncols=size(rho_sub,2);
                                current_pool=gcp('nocreate');
                                if isempty(current_pool)
                                    nsegments=1;
                                else
                                    nsegments=min(current_pool.NumWorkers,ncols);
                                end
                                edges=round(linspace(0,ncols,nsegments+1));
                                seg_starts=edges(1:end-1)+1;
                                seg_ends=edges(2:end);
                                fid_parts=cell(1,nsegments);

                                % Parallel processing over matrix segments
                                parfor seg=1:nsegments

                                    % Pull the current block of columns
                                    cols=seg_starts(seg):seg_ends(seg);
                                    rho_local=rho_sub(:,cols);
                                    cov_local=sparse(cols,1:numel(cols),1,ncols,numel(cols));
                                    fid_local=zeros([(nsteps+1) 1],'like',1i);

                                    % Loop over time steps
                                    for n=1:(nsteps+1)

                                        % Write local fid contribution
                                        fid_local(n)=hdot(coil_sub*cov_local,rho_local);

                                        % Step forward on rho and cov
                                        rho_local=P_sub*rho_local;
                                        cov_local=P_sub*cov_local;

                                    end

                                    % Store the segment contribution
                                    fid_parts{seg}=fid_local;

                                end

                                % Collect the results
                                for seg=1:nsegments
                                    answer=answer+fid_parts{seg};
                                end
                                
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

