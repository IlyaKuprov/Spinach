% Permutation symmetry treatment. Compiles character tables of 
% composite symmetry groups, builds the permutation table for
% each spin state in the basis, and builds projectors into the
% irreducible representations of the direct product symmetry 
% group. Syntax:
%
%             spin_system=symmetry(spin_system)
%
% Parameters:
%
%    spin_system  -  Spinach spin system description object
%                    produced as described in the spin system
%                    and basis specification sections, of the
%                    of the online manual.
%
% Outputs: symmetry factorisation information is written, for
%          each substance into the following fields
%
%              spin_system.bas.sym_fact(s).irr_dimensions
%              spin_system.bas.sym_fact(s).irr_projectors
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

function spin_system=symmetry(spin_system)

% Check the disable switch
if ismember('symmetry',spin_system.sys.disable)
    
    % Issue a reminder to the user
    report(spin_system,'WARNING - symmetry factorisation disabled by the user.');
    
    % Write empty cells and bounce out
    spin_system.bas.sym_a1g_only=NaN;
    spin_system.bas.sym_group={};
    spin_system.bas.sym_spins={}; return;
    
end

% The default is no symmetry
if ~isfield(spin_system.bas,'sym_group'), spin_system.bas.sym_group={}; end
if ~isfield(spin_system.bas,'sym_spins'), spin_system.bas.sym_spins={}; end
    
% Default A1g irrep switch
if ~isfield(spin_system.bas,'sym_a1g_only')

    % Depends on the formalism
    if strcmp(spin_system.bas.formalism,'zeeman-hilb')||...
       strcmp(spin_system.bas.formalism,'zeeman-wavef')

        % Hilbert space needs all irreps
        spin_system.bas.sym_a1g_only=false();

    else

        % Liouville space needs A1g only
        spin_system.bas.sym_a1g_only=true();

    end

end

% Check consistency
grumble(spin_system);
    
% Report back to the user if necessary
if ~isempty(spin_system.bas.sym_group)
    summary(spin_system,'symmetry','Symmetry summary');
end

% Loop over chemical substances
for s=1:spin_system.bas.nsubst

    % Get spin list and count for current substance
    spins_in_subst=spin_system.chem.parts{s}(:); 
    nspins_in_subst=numel(spins_in_subst);

    % Find out which symmetry groups belong to this substance
    groups_in_subst=false(1,numel(spin_system.bas.sym_group));
    for g=1:numel(spin_system.bas.sym_group)
        if all(ismember(spin_system.bas.sym_spins{g}(:),...
                        spins_in_subst))
            groups_in_subst(g)=true();
        end
    end

    % Pull out symmetry groups of the current substance
    sym_group_subst=spin_system.bas.sym_group(groups_in_subst);
    sym_spins_subst=spin_system.bas.sym_spins(groups_in_subst);

    % Translate indices to refer to the current substance
    for n=1:numel(sym_spins_subst)
        sym_spins_subst{n}=find(ismember(spins_in_subst,...
                                         sym_spins_subst{n}));
    end

    % Skip if the list is empty
    if isempty(sym_group_subst)
        spin_system.bas.sym_fact(s).irr_dimensions=[];
        spin_system.bas.sym_fact(s).irr_projectors={}; continue;
    end

    % Print substance header if some symmetry exists
    report(spin_system,['SYMMETRY TREATMENT, substance ' int2str(s)]);

    % Compute group direct product if necessary
    if numel(sym_group_subst)>1
        
        % Lift constituent groups from the database
        ngroups=numel(sym_group_subst); groups=cell(1,ngroups);
        for n=1:ngroups
            groups{n}=perm_group(sym_group_subst{n});
        end
        
        % Compute direct product character table
        group.characters=1;
        for n=1:ngroups
            group.characters=kron(group.characters,groups{n}.characters);
        end
        group.irrep_dims=group.characters(:,1)';
        group.n_irreps=size(group.characters,1);
        report(spin_system,['    ' num2str(group.n_irreps) ' irreps in the group direct product.']);
        report(spin_system,['    dimensions of the irreps ' num2str(group.irrep_dims)]);
        
        % Compute direct product element list
        group.elements=groups{1}.elements; group.order=groups{1}.order;
        for n=2:ngroups
            group.elements=[kron(group.elements,ones(groups{n}.order,1)) ...
                            kron(ones(group.order,1),groups{n}.elements+size(group.elements,2))];
            group.order=group.order*groups{n}.order;
        end
        report(spin_system,['    ' num2str(size(group.elements,1)) ' elements in the group direct product.']);
        
    else
        
        % Lift the group from the database
        group=perm_group(sym_group_subst{1});
        report(spin_system,['    ' num2str(group.n_irreps) ' irreps in the symmetry group.']);
        report(spin_system,['    dimensions of the irreps ' num2str(group.irrep_dims)]);
    
    end

    % Concatenate spin lists
    spins=horzcat(sym_spins_subst{:});
    
    % Run the SALC procedure
    if exist('group','var')
        
        % Pull out the current substance basis
        substance_basis=spin_system.bas.basis{s};
        bas_dimension=size(substance_basis,1);

        % Decide permutation table data type
        perm_tab_type=min_int_type(bas_dimension,'unsigned');

        % Preallocate basis set permutation table
        permutation_table=zeros(bas_dimension,group.order,perm_tab_type);

        % Compute the permutation table
        report(spin_system,'    computing the permutation table...');
        parfor n=1:group.order %#ok<*PFBNS>

            % Sequentially numbered spins
            group_element=1:nspins_in_subst;

            % Permutation of symmetry-related spins by a group element
            group_element(spins)=group_element(spins(group.elements(n,:)));

            % Corresponding permutation of basis set columns
            permuted_basis=substance_basis(:,group_element);

            % Find where symmetry-transformed states are in the old basis
            [~,index]=sortrows(permuted_basis); permuted_basis=[]; %#ok<NASGU>
            permutation_table(:,n)=cast(index,perm_tab_type);

        end

        % Clear a large array
        clear('substance_basis');

        % Compute irreducible representations
        if spin_system.bas.sym_a1g_only
            
            % Tell the user that only one irrep is pertinent
            report(spin_system,'    Liouville space mode - A1g only.');
            report(spin_system,'    building the fully symmetric irrep...');
            
            % Sort and prune the permutation table
            if (~isworkernode)&&(bas_dimension>1e4)
                
                % Distributed sorting of rows
                permutation_table=distrib_dim(permutation_table,1);
                permutation_table=sort(permutation_table,2,'ascend');
                permutation_table=gather(permutation_table);

            else

                % Single-thread sorting of rows
                permutation_table=sort(permutation_table,2,'ascend');

            end

            % Not sparse, so unihash() is not needed
            permutation_table=unique(permutation_table,'rows');

            % Get symmetrised subspace dimension
            irr_dimension=size(permutation_table,1);

            % Decide sparse constructor index data type
            idx_data_type=min_int_type(irr_dimension,'unsigned');

            % Construct the projector into the symmetrised subspace
            state_list=cast(1:irr_dimension,idx_data_type); state_list=state_list(:);
            index=[repmat(state_list,group.order,1) permutation_table(:)]; 
            clear('permutation_table'); index=unique(index,'rows');
            irr_projector=sparse(index(:,1),index(:,2),... % this needs
                                 ones(size(index,1),1),... % to be FP64
                                 irr_dimension,bas_dimension); clear('index');

            % Normalise the rows of the projector
            row_norms=sqrt(sum(irr_projector.^2,2));
            irr_projector=irr_projector./row_norms;

            % Convert into a single cell
            irr_projector={irr_projector};
            
            % Report the symmetrised subspace dimension back to the user
            report(spin_system,['    A1g irrep, ' int2str(irr_dimension) ' states.']);
            
        else
            
            % Tell the user that all irreps will be included
            report(spin_system,'    Hilbert space mode - all irreps included.');
            report(spin_system,'    building irreducible representations...');

            % Preallocate pertinent arrays
            irr_dimension=nan(group.n_irreps,1);
            irr_projector=cell(group.n_irreps,1);
            
            % Loop over irreps
            for n=1:group.n_irreps

                % Make SALC coefficient matrix element index
                vals=repmat(group.characters(n,:),bas_dimension,1);
                cols=cast(1:bas_dimension,perm_tab_type);
                cols=repmat(cols(:),1,group.order);
                rows=permutation_table;

                % Prune SALC coeff matrix element index
                rows=rows(:); cols=cols(:); vals=vals(:);
                kill_list=(vals==0); rows(kill_list)=[];
                cols(kill_list)=[];  vals(kill_list)=[];

                % Build SALC coefficient matrix
                coeff_matrix=sparse(rows,cols,vals);
                clear('rows','cols','vals');

                % Keep only non-zero SALCs in the projector
                coeff_matrix=coeff_matrix(:,any(coeff_matrix,1));
                
                % For non-empty irreps
                if nnz(coeff_matrix)>0

                    % Row index of the first NZ in each column
                    [~,fnz_row]=max(spones(coeff_matrix),[],1);
                    
                    % Convert to linear index of the first non-zero in each column
                    fnz_lin=sub2ind(size(coeff_matrix),fnz_row,1:size(coeff_matrix,2));
                    
                    % Sign of first NZ in each column
                    signs=sign(coeff_matrix(fnz_lin));

                    % Remove SALC sign ambiguity
                    coeff_matrix=signs.*coeff_matrix;

                    % Remove identical SALCs (columns)
                    coeff_matrix=transpose(coeff_matrix);
                    coeff_matrix=unihash(coeff_matrix);
                    coeff_matrix=transpose(coeff_matrix);

                end
                
                % Decide how to proceed
                if size(coeff_matrix,2)==0
                    
                    % Do nothing, this irrep is empty
                    
                % SALC orthogonalisation
                elseif group.irrep_dims(n)>1
                    
                    % Get the overlap connectivity matrix
                    overlap=coeff_matrix'*coeff_matrix;
                    overlap=(abs(overlap)>spin_system.tols.liouv_zero);
                    
                    % Find the target subspaces
                    member_states=scomponents(overlap); 
                    n_subspaces=max(member_states);
                    
                    % Preallocate the result
                    orth_coeff_matrix=cell(1,n_subspaces);
                    
                    % Fill the result
                    for k=1:n_subspaces

                        % Get subspace index
                        idx=(member_states==k);

                        % Inspect the subspace
                        if nnz(idx)==1

                            % Normalise one-vector SALCs
                            salc_vector=coeff_matrix(:,idx);
                            orth_coeff_matrix{k}=salc_vector/norm(salc_vector,2);

                        else
                            
                            % Orthonormalise bigger SALCs
                            vectors=full(coeff_matrix(:,idx));
                            vectors=clean_up(spin_system,orth(vectors),...
                                             spin_system.tols.liouv_zero);
                            orth_coeff_matrix{k}=sparse(vectors);

                        end

                    end
                    
                    % Reassemble the SALC coefficient matrix
                    coeff_matrix=cell2mat(orth_coeff_matrix);
                    
                else
                    
                    % For one-dimensional irreps, just normalize all SALCs
                    coeff_matrix=coeff_matrix./sqrt(sum(coeff_matrix.^2,1));
                    
                end
                
                % Assign projector arrays
                irr_projector{n}=coeff_matrix; 
                irr_dimension(n)=size(coeff_matrix,2);

                % Inform the user 
                report(spin_system,['    irrep #' int2str(n) ', ' ...
                                    int2str(group.irrep_dims(n)) '-dimensional, ' ...
                                    int2str(irr_dimension(n))    ' states.']);
                
            end
            
            % Remove zero-dimensional irreps
            keep_mask=logical(irr_dimension); 
            irr_dimension(~keep_mask)=[];
            irr_projector(~keep_mask)=[];

            % Let the user know
            if nnz(~keep_mask)>0
                report(spin_system,'    zero-dimensional irrep(s) removed.');
            end
            
        end
        
    end

    % Store symetry factorisations in the data structure
    spin_system.bas.sym_fact(s).irr_dimensions=irr_dimension;
    spin_system.bas.sym_fact(s).irr_projectors=irr_projector;
    
end

end

% Consistency enforcement
function grumble(spin_system)

% Check data types
if ~iscell(spin_system.bas.sym_group)
    error('bas.sym_group must be a cell array of strings.');
end
if ~iscell(spin_system.bas.sym_spins)
    error('bas.sym_spins must be a cell array of vectors.');
end

% Check dimensions
if numel(spin_system.bas.sym_spins)~=numel(spin_system.bas.sym_group)
    error('bas.sym_group and bas.sym_spins must have the same number of elements.');
end

% For non-empty specifications
if ~isempty(spin_system.bas.sym_group)

    % Check contents
    if any(~cellfun(@ischar,spin_system.bas.sym_group),'all')
        error('bas.sym_group must be a cell array of strings.');
    end
    if any(~cellfun(@isnumeric,spin_system.bas.sym_spins),'all')
        error('bas.sym_spins must be a cell array of vectors.');
    end

    % Check the spin indices
    for m=1:numel(spin_system.bas.sym_spins)

        % Check for common sense 
        if any(spin_system.bas.sym_spins{m}>spin_system.comp.nspins,'all')||...
           any(spin_system.bas.sym_spins{m}<1,'all')||...
           any(mod(spin_system.bas.sym_spins{m},1)~=0,'all')||...
              (numel(spin_system.bas.sym_spins{m})<2)
            error('incorrect spin labels in bas.sym_spins.');
        end

        % Check for group boundary violations
        for n=1:numel(spin_system.bas.sym_spins)
            if (n~=m)&&(~isempty(intersect(spin_system.bas.sym_spins{m},...
                                           spin_system.bas.sym_spins{n})))
                error('same spin is listed in multiple symmetry groups.');
            end
        end

        % Check for chemical boundary violations
        for n=1:numel(spin_system.chem.parts)

            % Get the intersection of symmetry and chemistry 
            common_spins=intersect(spin_system.bas.sym_spins{m}(:),...
                                   spin_system.chem.parts{n}(:));
            
            % Either all spins should be in a given substance, or none
            if (numel(common_spins)~=numel(spin_system.bas.sym_spins{m}))&&...
               (numel(common_spins)~=0)
                error('symmetries violate chemical partitions.');
            end

        end

    end
    
    % Check the group names
    for n=1:numel(spin_system.bas.sym_group)
        if ~ismember(spin_system.bas.sym_group{n},...
                     {'S2','S3','S4','S4A','S5','S6'})
            error('a group specified in bas.sym_group is not available.');
        end
    end

end

% Check the A1g irrep switch
if (~isnumeric(spin_system.bas.sym_a1g_only))&&...
   (~islogical(spin_system.bas.sym_a1g_only))&&...
   ((spin_system.bas.sym_a1g_only~=1)||...
    (spin_system.bas.sym_a1g_only~=0))
    error('the allowed values for bas.sym_a1g_only are 0 and 1.');
end

end

% I am regularly asked what an average Internet user can
% do to ensure his security. My first answer is usually
% "nothing, you're screwed".
%
% Bruce Schneier

