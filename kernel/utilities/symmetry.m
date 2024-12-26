% Permutation symmetry treatment. Compiles character tables of 
% composite symmetry groups, builds the permutation table for
% each spin state in the basis, and builds projectors into the
% irreducible representations of the direct product symmetry 
% group. Syntax:
%
%             spin_system=symmetry(spin_system,bas)
%
% Parameters:
%
%    spin_system  -  Spinach spin system description object
%                    produced as described in the spin system
%                    and basis specification sections, of the
%                    of the online manual.
%
%    bas          -  basis input structure described in the
%                    basis specification section of the manual
%
% Outputs:
%
%    spin_system.bas.irrep(n).projector - projector matrices
%                                         into each irreducible
%                                         representation
%
%    spin_system.bas.irrep(n).dimension - dimension of each ir-
%                                         reducible representa-
%                                         tion
%
% Note: this is a service function of the Spinach kernel that
%       should not be called directly; it is called by basis.m
%
% Note: non-Abelian groups and multi-dimensional irreps are sup-
%       ported - edit perm_group.m to add your own groups.
%
% ilya.kuprov@weizmann.ac.il
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=symmetry.m>

function spin_system=symmetry(spin_system,bas)

% Check consistency
grumble(spin_system,bas);

% Check the disable switch
if ismember('symmetry',spin_system.sys.disable)
    
    % Issue a reminder to the user
    report(spin_system,'WARNING - symmetry factorization disabled by the user.');
    
    % Write empty cells
    spin_system.comp.sym_group={};
    spin_system.comp.sym_spins={};
    spin_system.comp.sym_a1g_only=true();
    
else

    % Symmetry group
    if isfield(bas,'sym_group')
        spin_system.comp.sym_group=bas.sym_group;
    else
        spin_system.comp.sym_group={};
    end
    
    % Symmetry-related spins
    if isfield(bas,'sym_spins')
        spin_system.comp.sym_spins=bas.sym_spins;
    else
        spin_system.comp.sym_spins={};
    end
    
    % Irreducible representation composition
    if isfield(bas,'sym_a1g_only')
        spin_system.comp.sym_a1g_only=bas.sym_a1g_only;
    elseif strcmp(spin_system.bas.formalism,'zeeman-hilb')
        spin_system.comp.sym_a1g_only=false();
    else
        spin_system.comp.sym_a1g_only=true();
    end
    
    % Report back to the user
    if ~isempty(spin_system.comp.sym_group)
        summary(spin_system,'symmetry','permutation symmetry summary');
    end
    
    % Compute group direct product if necessary
    if numel(spin_system.comp.sym_group)>1
        
        % Lift constituent groups from the database
        ngroups=numel(spin_system.comp.sym_group); groups=cell(1,ngroups);
        for n=1:ngroups
            groups{n}=perm_group(spin_system.comp.sym_group{n});
        end
        
        % Compute direct product character table
        group.characters=1;
        for n=1:ngroups
            group.characters=kron(group.characters,groups{n}.characters);
        end
        group.irrep_dims=group.characters(:,1)';
        group.n_irreps=size(group.characters,1);
        report(spin_system,['' num2str(group.n_irreps) ' irreps in the group direct product.']);
        report(spin_system,['dimensions of the irreps ' num2str(group.irrep_dims)]);
        
        % Compute direct product element list
        group.elements=groups{1}.elements; group.order=groups{1}.order;
        for n=2:ngroups
            group.elements=[kron(group.elements,ones(groups{n}.order,1))...
                kron(ones(group.order,1),groups{n}.elements+size(group.elements,2))];
            group.order=group.order*groups{n}.order;
        end
        report(spin_system,['' num2str(size(group.elements,1)) ' symmetry operations in the group direct product.']);
        
        % Concatenate spin lists
        spins=horzcat(spin_system.comp.sym_spins{:});
        
    elseif isscalar(spin_system.comp.sym_group)
        
        % Lift the group from the database
        spins=spin_system.comp.sym_spins{1};
        group=perm_group(spin_system.comp.sym_group{1});
        report(spin_system,['' num2str(group.n_irreps) ' irreps in the symmetry group.']);
        report(spin_system,['dimensions of the irreps ' num2str(group.irrep_dims)]);
        
    else
        
        % Remind the user that symmetry is not operational
        report(spin_system,'no symmetry information available.');
        
    end
    
    % Run the SALC procedure
    if exist('group','var')
        
        % Preallocate the permutation table
        permutation_table=zeros(size(spin_system.bas.basis,1),group.order);
        
        % Compute the permutation table
        parfor n=1:group.order %#ok<*PFBNS>
            group_element=1:spin_system.comp.nspins;
            group_element(spins)=group_element(spins(group.elements(n,:)));
            permuted_basis=spin_system.bas.basis(:,group_element);
            [~,index]=sortrows(permuted_basis);
            permutation_table(:,n)=index;
        end
        
        % Compute irreducible representations
        if spin_system.comp.sym_a1g_only
            
            % Inform the user
            report(spin_system,'Liouville space mode - fully symmetric irrep only.');
            
            % Prune the permutation table
            symmetry_related_states=unique(sort(permutation_table,2,'ascend'),'rows');
            dimension=size(symmetry_related_states,1);
            
            % Populate the coefficient matrix
            index=unique([kron(ones(group.order,1),(1:dimension)') symmetry_related_states(:) ones(dimension*group.order,1)],'rows');
            coeff_matrix=sparse(index(:,1),index(:,2),index(:,3),dimension,size(spin_system.bas.basis,1));
            
            % Normalize the coefficient matrix
            norms=sqrt(sum(coeff_matrix.^2,2));
            coeff_matrix=spdiags(norms.^(-1),0,dimension,dimension)*coeff_matrix;
            
            % Report back to the user
            report(spin_system,['A1g irrep, ' num2str(dimension) ' states.']);
            
            % Return the projector and dimension
            spin_system.bas.irrep.projector=coeff_matrix';
            spin_system.bas.irrep.dimension=dimension;
            
        else
            
            % Inform the user
            report(spin_system,'full symmetry treatment - all irreps will be included.');
            report(spin_system,'processing irreducible representations...');
            
            % Determine the problem dimension
            basis_dim=size(spin_system.bas.basis,1);
            
            % Loop over irreducible representations
            for n=1:group.n_irreps
                
                % Build the transformation matrix
                rows=permutation_table;
                cols=kron((1:basis_dim)',ones(1,group.order));
                vals=kron(ones(basis_dim,1),group.characters(n,:));
                coeff_matrix=sparse(rows(:),cols(:),vals(:));
                clear('rows','cols','vals');
                
                % Remove sign ambiguity and clean up
                for k=1:basis_dim
                    non_zero_elements=nonzeros(coeff_matrix(:,k));
                    if any(non_zero_elements)
                        coeff_matrix(:,k)=coeff_matrix(:,k)*sign(non_zero_elements(1)); %#ok<SPRIX>
                    end
                end
                coeff_matrix=clean_up(spin_system,coeff_matrix,spin_system.tols.liouv_zero);
                
                % Remove zero columns
                hit_index=(sum(abs(coeff_matrix),1)==0); coeff_matrix(:,hit_index)=[];
                
                % Remove identical columns
                coeff_matrix=unique(coeff_matrix','rows')';
                
                % Decide how to proceed
                if size(coeff_matrix,2)==0
                    
                    % Do nothing
                    
                elseif group.irrep_dims(n)>1
                    
                    % Get the overlap matrix
                    overlap=logical(coeff_matrix'*coeff_matrix);
                    
                    % Find non-orthogonal subspaces
                    member_states=scomponents(overlap); n_subspaces=max(member_states);
                    
                    % Preallocate the result
                    orth_coeff_matrix=cell(1,n_subspaces);
                    
                    % Fill the result
                    for k=1:n_subspaces
                        vectors=full(coeff_matrix(:,member_states==k));
                        vectors=clean_up(spin_system,orth(vectors),spin_system.tols.liouv_zero);
                        orth_coeff_matrix{k}=sparse(vectors);
                    end
                    
                    % Build the matrix
                    coeff_matrix=cell2mat(orth_coeff_matrix);
                    
                else
                    
                    % Just normalize the SALCs
                    coeff_matrix=coeff_matrix./sqrt(sum(coeff_matrix.^2,1));
                    
                end
                
                % Inform the user and write the irrep into the data structure
                report(spin_system,['irreducible representation #' num2str(n) ...
                                    ', ' num2str(group.irrep_dims(n)) '-dimensional, ' ...
                                         num2str(size(coeff_matrix,2)) ' states.']);
                spin_system.bas.irrep(n).projector=coeff_matrix;
                spin_system.bas.irrep(n).dimension=size(coeff_matrix,2);
                
            end
            
            % Remove zero-dimensional irreps
            kill_mask=([spin_system.bas.irrep.dimension]==0);
            spin_system.bas.irrep(kill_mask)=[];
            if nnz(kill_mask)>0
                report(spin_system,'zero-dimensional irreps removed.');
            end
            
        end
        
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,bas)

% Check symmetry parameters
if isfield(bas,'sym_group')
    
    % Check the type
    if (~iscell(bas.sym_group))||any(~cellfun(@ischar,bas.sym_group))
        error('bas.sym_group must be a cell array of strings.');
    end
    
    % Check that bas.sym_spins exists
    if ~isfield(bas,'sym_spins')
        error('bas.sym_spins must be specified alongside bas.sym_group.');
    end
    
    % Check the type
    if (~iscell(bas.sym_spins))||any(~cellfun(@isnumeric,bas.sym_spins))
        error('bas.sym_spins must be a cell array of vectors.');
    end
    
    % Check the dimensions
    if numel(bas.sym_spins)~=numel(bas.sym_group)
        error('bas.sym_group and bas.sym_spins arrays must have the same number of elements.');
    end
    
    % Check the spin indices
    for m=1:length(bas.sym_spins)
        
        % Check for sense
        if any(bas.sym_spins{m}>spin_system.comp.nspins)||any(bas.sym_spins{m}<1)||(numel(bas.sym_spins{m})<2)
            error('incorrect spin labels in bas.sym_spins.');
        end
        
        % Check for intersections
        for n=1:length(bas.sym_spins)
            if (n~=m)&&(~isempty(intersect(bas.sym_spins{m},bas.sym_spins{n})))
                error('same spin is listed in multiple symmetry groups in bas.sym_spins.');
            end
        end
        
    end
    
    % Check the group names
    for n=1:length(bas.sym_group)
        if ~ismember(bas.sym_group{n},{'S2','S3','S4','S4A','S5','S6'})
            error('the group requested in bas.sym_group is not available.');
        end
    end
    
    % Check the irrep switch
    if isfield(bas,'sym_a1g_only')&&(~isnumeric(bas.sym_a1g_only))&&(~islogical(bas.sym_a1g_only))&&((bas.sym_a1g_only~=1)||(bas.sym_a1g_only~=0))
        error('the allowed values for bas.sym_a1g_only are 0 and 1.');
    end
    
else
    
    % Enforce no sys.sym_spins without sys.sym_group
    if isfield(bas,'sym_spins')
        error('bas.sym_group must be specified alongside bas.sym_spins.');
    end
    
    % Enforce no sys.sym_a1g_only without sys.sym_group
    if isfield(bas,'sym_a1g_only')
        error('bas.sym_group must be specified alongside bas.sym_a1g_only.');
    end
    
end

end

% I am regularly asked what an average Internet user can
% do to ensure his security. My first answer is usually
% "nothing, you're screwed".
%
% Bruce Schneier

