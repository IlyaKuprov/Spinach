% Liouvillian path tracing. Treats the user-supplied Liouvillian 
% as the adjacency matrix of a graph, computes the weakly connect-
% ed subgraphs of that graph and returns a cell array of project-
% ors into independently evolving populated subspaces. Syntax:
%
%            projectors=path_trace(spin_system,L,rho)
%
% Parameters:
%
%     L   -  Hamiltonian or Liouvillian matrix 
%
%     rho -  the initial state (source state screening)
%            or the detection state (destination state
%            screening); pass [] to disable screening
%
% Outputs:
%
%   projectors - a cell array of projectors into independently
%                evolving populated subspaces. The projectors 
%                are to be used as follows:
%
%                   L_reduced=P'*L*P;    (for matrices)
%                   rho_reduced=P'*rho;  (for state vectors)
%
% Note: further information on how this function works is availa-
%       ble in our papers on this subject
%
%             http://dx.doi.org/10.1063/1.3398146
%             http://dx.doi.org/10.1016/j.jmr.2011.03.010
%
% i.kuprov@soton.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=path_trace.m>

function projectors=path_trace(spin_system,L,rho)

% Check the input
grumble(L,rho);

% Check run conditions
if ismember('pt',spin_system.sys.disable)
    
    % Return a unit projector if path tracing is disabled
    report(spin_system,'WARNING - path tracing disabled by the user.');
    projectors={1}; return
    
elseif size(L,2)<spin_system.tols.merge_dim
    
    % Return a unit projector if the space is small anyway
    report(spin_system,'small state space - path tracing skipped.');
    projectors={1}; return
    
else
    
    % Report to the user
    report(spin_system,['analyzing ' num2str(size(L,1)) '-dimensional state space.']);
    report(spin_system,['Liouvillian zero tolerance ' num2str(spin_system.tols.liouv_zero)]);
    report(spin_system,['population tolerance for subspace drop ' num2str(spin_system.tols.subs_drop)]);
    
end

% Get the connectivity matrix
G=(abs(L)>spin_system.tols.liouv_zero);

% Make sure isolated states do not get lost
G=or(G,transpose(G)); G=or(G,speye(size(G)));

% Get the weakly connected subgraphs
member_states=scomponents(G);

% Determine the number of subspaces
n_subspaces=max(member_states);
report(spin_system,['found ' num2str(n_subspaces) ' non-interacting subspaces.']);
report(spin_system,'running subspace population analysis...');

% All subspaces matter by default
subspace_important=true(n_subspaces,1);

% Screen by population
if ~isempty(rho)
    
    % Extract switches for parfor loop
    tolerance=spin_system.tols.subs_drop;
    formalism=spin_system.bas.formalism;
    
    % Look at actual populations
    parfor n=1:n_subspaces
        
        % Determine the importance
        if ismember(formalism,{'sphten-liouv','zeeman-liouv'})
            
            % Liouville space has state vectors
            subspace_important(n)=(norm(rho.*(member_states==n),1)>tolerance);
            
        elseif ismember(formalism,{'zeeman-hilb'})
            
            % Hilbert space has state matrices
            subspace_important(n)=(norm(rho(member_states==n,:),1)>tolerance)|...
                                  (norm(rho(:,member_states==n),1)>tolerance);
                              
        else
            
            % Complain and bomb out
            error('unexpected formalism specification.');
            
        end
        
    end
    
end

% Ignore unpopulated subspaces
significant_subspaces=find(subspace_important);
n_subspaces=numel(significant_subspaces);

% Preallocate projectors and counters
projectors=cell(1,n_subspaces);

% Build projectors into significant subspaces
parfor n=1:n_subspaces
    
    % Find states populating the current subspace
    state_index=find(member_states==significant_subspaces(n));

    % Determine subspace dimension
    subspace_dim=numel(state_index);
    
    % Build the projector into the current subspace
    projectors{n}=sparse(state_index,1:subspace_dim,ones(1,subspace_dim),size(L,1),subspace_dim);
                     
end

% Compile dimension statistics
subspace_dims=cellfun(@(x)size(x,2),projectors);
unique_dims=unique(subspace_dims);

% Report to the user
for n=1:numel(unique_dims)
    current_dim_count=nnz(subspace_dims==unique_dims(n));
    report(spin_system,[num2str(current_dim_count)...
                        ' populated subspaces of dimension '...
                        num2str(unique_dims(n))]);
end
report(spin_system,['keeping a total of ' num2str(n_subspaces) ' independent subspaces '...
                    'of total dimension ' num2str(sum(cellfun(@(x)size(x,2),projectors)))]);
                
% Merge small subspaces
if ismember('merge',spin_system.sys.disable)
    
    % Inform the user
    report(spin_system,'WARNING - small subspace merging disabled by the user.');

else
    
    % Inform the user
    report(spin_system,['merging small subspaces into batches of dimension '...
                         num2str(spin_system.tols.merge_dim) '...']);
    
    % Call the bin packer
    bins=binpack(subspace_dims,spin_system.tols.merge_dim);
    
    % Group the subspaces
    new_projectors=cell(numel(bins),1);
    for n=1:numel(bins)
        new_projectors{n}=[projectors{bins{n}}];
    end
    projectors=new_projectors;
    
    % Report to the user
    for n=1:numel(projectors)
        report(spin_system,['working subspace ' num2str(n) ', dimension '  num2str(size(projectors{n},2))]);
    end
    
end

end

% Consistency enforcement
function grumble(L,rho)
if (~isnumeric(L))||(~isnumeric(rho))
    error('both inputs must be numeric.');
end
if size(L,1)~=size(L,2)
    error('L must be square.');
end
if (~isempty(rho))&&(size(L,2)~=size(rho,1))
    error('dimensions of L and rho must be consistent.');
end
end

% "My dear fellow, who will let you?"
% "That's not the point. The point is, who will stop me?"
%
% Ayn Rand, "The Fountainhead"

