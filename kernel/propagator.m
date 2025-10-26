% Calculates exponential propagator exp(-1i*L*timestep) using scaled 
% and squared Taylor series method. Syntax:
%
%                 P=propagator(spin_system,L,timestep)
%
% Parameters:
%
%    L          -  Hamiltonian or Liouvillian matrix
%
%    timestep   -  propagation time step
%
% Outputs:
%
%    P          -  propagator matrix
%
% Note: GPUs are supported, add 'gpu' to sys.enable array during 
%       calculation setup.
%
% Note: propagator caching (https://doi.org/10.1063/1.4928978) is
%       supported, add 'prop_cache' to sys.enable array to enable.
%
% Note: we did have Chebyshev and Newton series here at one point,
%       as well as the Pade method. None of them had lived up to
%       their marketing.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=propagator.m>

function P=propagator(spin_system,L,timestep)

% Check consistency
grumble(L,timestep); tic;

% Fast bypass for small matrices
if size(L,1)<spin_system.tols.small_matrix
    P=expm((-1i*timestep)*L); return; 
end

% Check the cache
if ismember('prop_cache',spin_system.sys.enable)

    % Hash the generator, the step, and the tolerance
    prop_hash=md5_hash({L,timestep,spin_system.tols.prop_chop});
    
    % Generate the cache record name in the scratch directory
    filename=[spin_system.sys.scratch filesep 'spinach_prop_' prop_hash '.mat'];
    
    % Check if the file exists
    if exist(filename,'file')

        % Try to use
        try
            
            % Try to load the cache record
            load(filename,'P'); load_time=toc();
            
            % Check load success
            if exist('P','var')
                
                % Report the stats and return the cached propagator
                report(spin_system,['cache record loaded in ' num2str(load_time) ' seconds']);  
                report(spin_system,['propagator dimension ' num2str(size(P,1)) ...
                                    ', nnz ' num2str(nnz(P))                   ...
                                    ', density ' num2str(100*nnz(P)/numel(P))  ...
                                    '%, sparsity ' num2str(issparse(P))]);
                return;

            else
                
                % Do not make a fuss on fail
                report(spin_system,'could not read the cache record, recomputing...');
                
            end
            
        catch

            % Do not make a fuss on fail
            report(spin_system,'could not read the cache record, recomputing...');

        end

    else
        
        % Let the user know the propagator is being computed
        report(spin_system,'cache record not found, computing...');

    end
    
end

% Set and clean up a shorthand for -i*L*dt
A=clean_up(spin_system,(-1i*timestep)*L,spin_system.tols.prop_chop);

% Inform the user about generator densities
report(spin_system,['generator dimension ' num2str(size(A,1)) ...
                    ', nnz ' num2str(nnz(A))                  ...
                    ', density ' num2str(100*nnz(A)/numel(A)) ...
                    '%, sparsity ' num2str(issparse(A))]);

% Estimate the norm
mat_norm=cheap_norm(A);

% Check the norm
if mat_norm>1e9
    
    % The user is doing something silly, bomb out
    error('norm of -i*L*dt exceeds 1e9, check your L and your dt.');

elseif mat_norm>1024
    
    % The user is really pushing it, take precautionary measures
    report(spin_system,'WARNING - time step greatly exceeds the timescale of system dynamics.');
    report(spin_system,'WARNING - exponentiation tolerance tightened up to eps(''double'').');
    spin_system.tols.prop_chop=eps('double');
    
elseif mat_norm>16
    
    % Inform the user just in case
    report(spin_system,'WARNING - time step exceeds the timescale of system dynamics.');
    
end

% Determine scaling and squaring parameters
n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings;
report(spin_system,['scaling -i*L*dt down by ' num2str(scaling_factor) ...
                    ' and squaring the propagator ' num2str(n_squarings) ' times.']);

% Scale the matrix
if scaling_factor>1, A=(1/scaling_factor)*A; end

% Get the propagator
if ismember('gpu',spin_system.sys.enable)&&(size(A,1)>500)
    
    % Run Taylor series procedure on the GPU
    A=gpuArray(A); P=speye(size(A));
    next_term=gpuArray.speye(size(A)); n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+gather(next_term); n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on GPU in ' num2str(n) ' iterations.']);
    
else
    
    % Run Taylor series procedure on the CPU
    P=speye(size(A)); next_term=P; n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+next_term; n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on CPU in ' num2str(n) ' iterations.']);
    
end

% Reclaim memory
clear('A','next_term');

% Clean up the result
P=clean_up(spin_system,P,spin_system.tols.prop_chop);

% Inform the user
report(spin_system,['propagator dimension ' num2str(size(P,1)) ...
                    ', nnz ' num2str(nnz(P))                   ...
                    ', density ' num2str(100*nnz(P)/numel(P))  ...
                    '%, sparsity ' num2str(issparse(P))]);

% Run the squaring stage
if n_squarings>0
    
    % Run the appropriate squaring process
    if ismember('gpu',spin_system.sys.enable)
        
        % Move to GPU
        P=gpuArray(P);
        
        % Run squaring on GPU
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['GPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
        
        % Gather the propagator
        P=gather(P);
        
    else
        
        % Run CPU squaring
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['CPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
      
    end
    
    % Report matrix dimension and density statistics
    report(spin_system,['propagator dimension ' num2str(size(P,1)) ...
                        ', nnz ' num2str(nnz(P))                   ...
                        ', density ' num2str(100*nnz(P)/numel(P))  ...
                        '%, sparsity ' num2str(issparse(P))]);

end
    
% Write the cache record if caching is beneficial
if ismember('prop_cache',spin_system.sys.enable)&&(toc>0.01)

    % Do not fight other workers
    if ~exist(filename,'file')
    
        % Try to save
        try
            
            % Modern format, compressed
            save(filename,'P','-v7.3'); drawnow;
            report(spin_system,'propagator cache record saved.');

        catch

            % Do not make a fuss on fail, this can happen
            % for large parallel pools where many workers
            % may be trying to write the same file.
            report(spin_system,'could not save cache record.');

        end

    end
    
elseif ismember('prop_cache',spin_system.sys.enable)
    
    % Tell the user that caching is pointless here
    report(spin_system,'cache record not worth saving.');
    
end

end

% Consistency enforcement
function grumble(L,timestep)
if isa(L,'polyadic')
    error('L is a polyadic - use evolution() instead.');
end
if (~isnumeric(L))||(size(L,1)~=size(L,2))
    error('L argument must be a square matrix.');
end
if (~isnumeric(timestep))||(~isscalar(timestep))||...
   (~isfinite(timestep))
    error('timestep must be a finite scalar.');
end
end

% To preserve one's mind intact through a modern college education is a
% test of courage and endurance, but the battle is worth it and the sta-
% kes are the highest possible to man: the survival of reason.
%
% Ayn Rand

