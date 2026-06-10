% Basis set generation. This is the second mandatory function 
% (after create.m) that must be called in every calculation to
% build the spin_system data structure. Syntax:
%
%              spin_system=basis(spin_system,bas)
%
% Parameters:
%
%     spin_system  - Spinach data structure, the output of
%                    create.m function
%
%     bas          - basis set specification structure de-
%                    scribed in the online manual
%
% Outputs:
%
%     spin_system  - Spinach data structure, updated with
%                    the basis set information
%
% Note: it is important to understand the basis set selection 
%       process in spin dynamics simulations, see
%
%           https://doi.org/10.1007/978-3-031-05607-9_7
%           https://doi.org/10.1063/1.3624564
%           https://doi.org/10.1002/mrc.4660
%
%       for further information. For data structure layouts,
%       see kernel/conventions/object_diagrams folder.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=basis.m>

function spin_system=basis(spin_system,bas)

% Show the basis module banner
banner(spin_system,'basis_banner');

% Check consistency
grumble(spin_system,bas);

% Store settings in the structure
spin_system.bas=bas; clear('bas');

% Report back to the user
summary(spin_system,'basis_settings');

% Record the number of substances here too
spin_system.bas.nsubst=numel(spin_system.chem.parts);

% Preallocate basis descriptor array, per substance
spin_system.bas.basis=cell(spin_system.bas.nsubst,1);

% Preallocate spherical tensor characteristics
if strcmp(spin_system.bas.formalism,'sphten-liouv')

    % Total projection quantum number of each basis state
    spin_system.bas.tot_proj=cell(spin_system.bas.nsubst,1);   

    % Number of spins correlated by each basis state
    spin_system.bas.tot_cord=cell(spin_system.bas.nsubst,1);

end

% Decide integer type for spin state indexing
max_state_index=max(spin_system.comp.mults.^2-1);
state_idx_type=min_int_type(max_state_index,'signed');

% Loop over chemical substances
for s=1:spin_system.bas.nsubst

    % Get spin list and count for current substance
    spins_in_subst=spin_system.chem.parts{s}(:); 
    nspins_in_subst=numel(spins_in_subst);

    % Isolate pertinent sub-arrays for the current substance
    mults_in_subst=spin_system.comp.mults(spins_in_subst);
    isots_in_subst=spin_system.comp.isotopes(spins_in_subst);
    inter_in_subst=spin_system.inter.coupling.matrix(spins_in_subst,...
                                                     spins_in_subst);

    % Proximity matrix had been computed in create.m
    prxmat_in_subst=spin_system.inter.proxmatrix(spins_in_subst,...
                                                 spins_in_subst);

    % Process spherical tensor basis sets
    if strcmp(spin_system.bas.formalism,'sphten-liouv')

        % Connectivity analysis for IK-DNP basis set
        if strcmp(spin_system.bas.approximation,'IK-DNP')

            % Find electrons and nuclei in the current substance
            e_idx=cellfun(@iselectron,isots_in_subst);
            n_idx=cellfun(@isnucleus, isots_in_subst);

            % Make sure there are only electrons and nuclei
            if (nnz(e_idx)==0)||(nnz(n_idx)==0)
                error('IK-DNP approximation requires both electrons and nuclei.');
            end
            if ~all(e_idx|n_idx,'all')
                error('IK-DNP approximation can only handle electrons and nuclei.');
            end

            % Isolate three distinct types of interactions (inter-electron, electron-nucl, inter-nucl)
            ee_couplings=inter_in_subst; ee_couplings(:    ,n_idx)={[]}; ee_couplings(n_idx,:    )={[]};
            en_couplings=inter_in_subst; en_couplings(e_idx,e_idx)={[]}; en_couplings(n_idx,n_idx)={[]};
            nn_couplings=inter_in_subst; nn_couplings(:    ,e_idx)={[]}; nn_couplings(e_idx,:    )={[]};

            % Generate three types of connectivity matrices (inter-electron, electron-nucl, inter-nucl)
            ee_conmat=sparse(cellfun(@(x)norm(x,2),ee_couplings)>2*pi*spin_system.tols.inter_cutoff);
            en_conmat=sparse(cellfun(@(x)norm(x,2),en_couplings)>2*pi*spin_system.tols.inter_cutoff);
            nn_conmat=sparse(cellfun(@(x)norm(x,2),nn_couplings)>2*pi*spin_system.tols.inter_cutoff);

            % Clear large temporary arrays
            clear('inter_in_subst','ee_couplings',...
                  'en_couplings','nn_couplings');

            % Make sure each spin is connected to itself
            ee_conmat=ee_conmat|speye(size(ee_conmat));
            en_conmat=en_conmat|speye(size(en_conmat));
            nn_conmat=nn_conmat|speye(size(nn_conmat));
  
            % Make sure connectivity is reciprocal
            ee_conmat=ee_conmat|transpose(ee_conmat);
            en_conmat=en_conmat|transpose(en_conmat);
            nn_conmat=nn_conmat|transpose(nn_conmat);

            % Issue a report to the user on the current substance
            report(spin_system,['chemical substance ' int2str(s) ':']);
            report(spin_system,['    EE connectivity matrix density: ' ...
                                num2str(100*nnz(ee_conmat)/numel(ee_conmat)) '%']);
            report(spin_system,['    EN connectivity matrix density: ' ...
                                num2str(100*nnz(en_conmat)/numel(en_conmat)) '%']);
            report(spin_system,['    NN connectivity matrix density: ' ...
                                num2str(100*nnz(nn_conmat)/numel(nn_conmat)) '%']);
           
        end

        % Connectivity analysis for IK-1,2 basis sets
        if ismember(spin_system.bas.approximation,{'IK-1','IK-2'})
        
            % Build connectivity matrix
            switch spin_system.bas.connectivity
            
                case 'scalar_couplings'

                    % Use scalar parts of all interaction tensors
                    conmat_in_subst=sparse(abs(cellfun(@trace,inter_in_subst)/3)>...
                                           2*pi*spin_system.tols.inter_cutoff);
                
                case 'full_tensors'
                
                    % Use complete interaction tensors
                    conmat_in_subst=sparse(cellfun(@(x)norm(x,2),inter_in_subst)>...
                                           2*pi*spin_system.tols.inter_cutoff);

                otherwise

                    % Complain and bomb out
                    error('unrecognised connectivity type.');
                
            end

            % Clear a large array
            clear('inter_in_subst');
        
            % Make sure each spin is connected and proximate to itself
            conmat_in_subst=conmat_in_subst|speye(size(conmat_in_subst));
            prxmat_in_subst=prxmat_in_subst|speye(size(prxmat_in_subst));
        
            % Make sure connectivity and proximity are reciprocal
            conmat_in_subst=conmat_in_subst|transpose(conmat_in_subst);
            prxmat_in_subst=prxmat_in_subst|transpose(prxmat_in_subst);

            % Count disconnected subnetworks
            n_sub_con=max(scomponents(conmat_in_subst));
            n_sub_prx=max(scomponents(prxmat_in_subst));
            n_sub_all=max(scomponents(conmat_in_subst|prxmat_in_subst));
        
            % Issue a report to the user on the current substance
            report(spin_system,['CHEMICAL SUBSTANCE ' int2str(s) ':']);
            report(spin_system,['    connectivity matrix density: ' ...
                                num2str(100*nnz(conmat_in_subst)/numel(conmat_in_subst)) '%']);
            report(spin_system,['    number of partitions in CM:  ' int2str(n_sub_con)]);
            report(spin_system,['    proximity matrix density:    ' ...
                                num2str(100*nnz(prxmat_in_subst)/numel(prxmat_in_subst)) '%']);
            report(spin_system,['    number of partitions in PM:  ' int2str(n_sub_prx)]);
            report(spin_system,['    non-inter. subsystem count:  ' int2str(n_sub_all)]);
            
        end
    
        % Build state lists for individual spins
        spin_state_lists_in_subst=cell(nspins_in_subst,1);
        for n=1:nspins_in_subst

            % Spinach spin state indexing convention 
            state_list=0:(mults_in_subst(n)^2-1);
            state_list=cast(state_list,state_idx_type);
            spin_state_lists_in_subst{n}=state_list(:);

        end
    
        % Apply longitudinal state filters
        if isfield(spin_system.bas,'longitudinal')&&...
           (~isempty(spin_system.bas.longitudinal{s}))
            
            % No longitudinal state filters by default
            longit_spins_in_subst=false(1,nspins_in_subst);

            % Over specifications
            for k=1:numel(spin_system.bas.longitudinal{s})

                % Specification by spin number
                if isnumeric(spin_system.bas.longitudinal{s}{k})
                    
                    % Find the spin on the index of the current substance
                    which_spins=ismember(spins_in_subst,spin_system.bas.longitudinal{s}{k});
                    if any(which_spins,'all')

                        % Add to the longitudinal filter index
                        longit_spins_in_subst=longit_spins_in_subst|which_spins;

                        % Inform the user
                        report(spin_system,['    T(L,M=0) states only for spin ' ...
                                                 int2str(spin_system.bas.longitudinal{s}{k})]);

                    else

                        % Complain and bomb out
                        error(['spin ' int2str(spin_system.bas.longitudinal{s}{k}) ...
                               ' does not belong to substance ' int2str(s)]);

                    end
                
                % Specification by isotope string
                elseif ischar(spin_system.bas.longitudinal{s}{k})

                    % Find spins on the isotope list of the current substance
                    which_spins=ismember(isots_in_subst,spin_system.bas.longitudinal{s}(k));
                    if any(which_spins,'all')

                        % Add to the longitudinal filter index
                        longit_spins_in_subst=longit_spins_in_subst|which_spins;

                        % Inform the user
                        report(spin_system,['    T(L,M=0) states only for ' ...
                                                 spin_system.bas.longitudinal{s}{k}]);
                    else

                        % Complain and bomb out
                        error(['no instances of ' spin_system.bas.longitudinal{s}{k} ...
                               ' found in substance ' int2str(s)]);

                    end

                else

                    % Complain and bomb out
                    error('unrecognised longitudinal state specification.');

                end

            end

            % Loop over spins to be filtered
            for n=find(longit_spins_in_subst)

                % Remove T(L,M~=0) states from state lists
                [~,M]=lin2lm(spin_state_lists_in_subst{n}); 
                spin_state_lists_in_subst{n}(logical(M))=[];

            end
        
        end
    
        % Compute single-spin subspace dimensions after filters
        spin_dims_in_subst=cellfun(@numel,spin_state_lists_in_subst);

        % Generate subgraphs
        switch spin_system.bas.approximation
        
            case 'none'
            
                % Single subgraph with all spins in it
                cpl_subgraphs_in_subst=true(1,nspins_in_subst);
            
                % Do not run proximity analysis
                prx_subgraphs_in_subst=[];
            
            case 'IK-0'

                % Make sure basis set level on the coupling network makes sense
                actual_inter_level=min([nspins_in_subst spin_system.bas.inter_level]);
            
                % All possible sets of actual_level spins
                col_index=nchoosek(spins_in_subst,actual_inter_level);
                
                % Get the number of sets
                nsets=size(col_index,1);
                
                % Assign numbers to sets
                row_index=repmat((1:nsets)',1,actual_inter_level);
                
                % Generate combinatorial subgraph list
                cpl_subgraphs_in_subst=sparse(row_index,col_index, ...
                                              true,nsets,nspins_in_subst);

                % Skip proximity analysis
                prx_subgraphs_in_subst=[];

                % Inform the user
                report(spin_system,['    IK-0: ' int2str(nsets) ...
                                    ' combinatorial subgraphs']);
        
            case 'IK-1'

                % Make sure basis set levels make sense
                actual_inter_level=min([nspins_in_subst spin_system.bas.inter_level]);
                actual_prox_level=min([nspins_in_subst spin_system.bas.prox_level]);
        
                % Run coupling connectivity analysis and get subgraphs
                cpl_subgraphs_in_subst=dfpt(conmat_in_subst,actual_inter_level);

                % Inform the user
                report(spin_system,['    IK-1: ' int2str(size(cpl_subgraphs_in_subst,1)) ...
                                    ' coupling subgraphs']);
            
                % Run spatial proximity analysis and get subgraphs
                prx_subgraphs_in_subst=dfpt(prxmat_in_subst,actual_prox_level);

                % Inform the user
                report(spin_system,['    IK-1: ' int2str(size(prx_subgraphs_in_subst,1)) ...
                                    ' proximity subgraphs']);
        
            case 'IK-2'

                % Make sure basis set space level makes sense
                actual_prox_level=min([nspins_in_subst spin_system.bas.prox_level]);
        
                % Coupling subgraphs use connectivity 
                cpl_subgraphs_in_subst=unique(conmat_in_subst,'rows');
            
                % Inform the user
                report(spin_system,['    IK-2: ' int2str(size(cpl_subgraphs_in_subst,1)) ...
                                    ' coupling subgraphs']);
            
                % Run proximity analysis and get subgraphs
                prx_subgraphs_in_subst=dfpt(prxmat_in_subst,actual_prox_level);

                % Inform the user
                report(spin_system,['    IK-2: ' int2str(size(prx_subgraphs_in_subst,1)) ...
                                    ' proximity subgraphs']);

            case 'IK-DNP'

                % Make sure basis set levels make sense
                actual_ee_level=min([nspins_in_subst spin_system.bas.inter_level(1)]);
                actual_en_level=min([nspins_in_subst spin_system.bas.inter_level(2)]);
                actual_nn_level=min([nspins_in_subst spin_system.bas.inter_level(3)]);
                
                % Inter-electron connectivity analysis
                ee_subgraphs=dfpt(ee_conmat,actual_ee_level);

                % Inform the user, do not mention trivial subgraphs
                report(spin_system,['    IK-DNP: ' int2str(size(ee_subgraphs,1)-nnz(n_idx)) ...
                                    ' inter-electron subgraphs']);
               
                % Electron-nuclear connectivity analysis
                en_subgraphs=dfpt(en_conmat,actual_en_level);

                % Inform the user
                report(spin_system,['    IK-DNP: ' int2str(size(en_subgraphs,1)) ...
                                    ' electron-nuclear subgraphs']);

                % Inter-nuclear connectivity analysis
                nn_subgraphs=dfpt(nn_conmat,actual_nn_level);

                % Inform the user, do not mention trivial subgraphs
                report(spin_system,['    IK-DNP: ' int2str(size(nn_subgraphs,1)-nnz(e_idx)) ...
                                    ' inter-nuclear subgraphs']);

                % Merge coupling subgraph lists
                cpl_subgraphs_in_subst=[ee_subgraphs; en_subgraphs; nn_subgraphs];

                % Do not run proximity analysis
                prx_subgraphs_in_subst=[];

            otherwise
            
                % Complain and bomb out
                error('unrecognised basis set.');
        
        end

        % Include user-specified subgraphs
        if (~isfield(spin_system.bas,'manual'))||...
             isempty(spin_system.bas.manual)

            % Empty list if the user had said nothing
            man_subgraphs_in_subst=false(0,nspins_in_subst);

        else

            % Get the manual subgraph array started
            man_subgraphs_in_subst=false(0,nspins_in_subst);

            % Loop over manual subgraph specs
            for n=1:size(spin_system.bas.manual,1)

                % Project the subgraph into the current substance
                current_subgraph=spin_system.bas.manual(n,spins_in_subst);

                % Check chemical boundaries
                if nnz(current_subgraph)==0

                    % Belongs to another substance, do nothing

                elseif nnz(current_subgraph)==nnz(spin_system.bas.manual(n,:))

                    % Belongs to this substance, add to the list
                    man_subgraphs_in_subst=[man_subgraphs_in_subst;
                                            current_subgraph]; %#ok<AGROW>

                else

                    % Complain and bomb out
                    error('manual subgraph spec crosses chemical boundaries.');

                end

            end

            % Inform the user
            report(spin_system,['    added ' int2str(size(man_subgraphs_in_subst,1)) ...
                                ' manual subgraphs']);

        end

        % Assemble the full subgraph list
        subgraphs_in_subst=[cpl_subgraphs_in_subst; ...
                            prx_subgraphs_in_subst; ...
                            man_subgraphs_in_subst];
        subgraphs_in_subst=logical(subgraphs_in_subst);
    
        % Ignore spin zero particles
        spin_zero=(mults_in_subst==1);
        if nnz(spin_zero)>0
            report(spin_system,['    ' int2str(nnz(spin_zero)) ...
                                ' particles have zero spin']);
            subgraphs_in_subst(:,spin_zero)=0;
        end
    
        % Remove empty subgraphs
        empty_ones=~logical(sum(subgraphs_in_subst,2)); 
        subgraphs_in_subst(empty_ones,:)=[];

        % Remove identical and completely enclosed subgraphs
        report(spin_system,'    removing redundant subgraphs...');
        subgraphs_in_subst=unique(subgraphs_in_subst,'rows');
        subgraphs_in_subst=prune_subgraphs(subgraphs_in_subst);
        
        % Report subgraph sizes
        sg_sizes_in_subst=sum(subgraphs_in_subst,2);
        for n=min(sg_sizes_in_subst):max(sg_sizes_in_subst)
            if nnz(sg_sizes_in_subst==n)>0
                report(spin_system,['    keeping ' num2str(nnz(sg_sizes_in_subst==n)) ...
                                    ' subgraphs with ' num2str(n) ' spins']);
            end
        end

        % Report projection quantum numbers
        if isfield(spin_system.bas,'projections')&&...
           (~isempty(spin_system.bas.projections))&&...
           (~isempty(spin_system.bas.projections{s}))
            report(spin_system,['    keeping states with total M=[' ...
                                num2str(spin_system.bas.projections{s}) ']']);
            projections_subst=spin_system.bas.projections{s};
        else
            projections_subst=[];
        end

        % Decide zero-quantum subset
        if isfield(spin_system.bas,'zero_quantum')&&...
           (~isempty(spin_system.bas.zero_quantum))&&...
           (~isempty(spin_system.bas.zero_quantum{s}))
            
            % Localise spec to the current substance
            zq_spins_in_subst=false(1,nspins_in_subst);

            % Loop over the specification
            for n=1:numel(spin_system.bas.zero_quantum{s})

                % For spin number specification
                if isnumeric(spin_system.bas.zero_quantum{s}{n})
                    
                    % Update the spin list for the zero-quantum restriction
                    zq_spins_in_subst=zq_spins_in_subst|...
                                      ismember(spins_in_subst,spin_system.bas.zero_quantum{s}{n});

                % For isotope string specification
                elseif ischar(spin_system.bas.zero_quantum{s}{n})

                    % Update the spin list for the zero-quantum restriction
                    zq_spins_in_subst=zq_spins_in_subst|...
                                      ismember(isots_in_subst,spin_system.bas.zero_quantum{s}(n));

                else

                    % Complain and bomb out
                    error('unrecognised zero-quantum restriction specification.');

                end

            end

            % Print a report
            if any(zq_spins_in_subst,'all')
                report(spin_system,['    keeping ZQ correlations of spins: ' ...
                                         int2str(find(zq_spins_in_subst))]);
            end

        else

           % Empty declaration needed for parfor later
           zq_spins_in_subst=false(1,nspins_in_subst); 

        end

        % Balance subgraph list for parallel processing
        subgraphs_in_subst=logical(subgraphs_in_subst);
        shuff=randperm(size(subgraphs_in_subst,1));
        subgraphs_in_subst=subgraphs_in_subst(shuff,:);

        % Populate basis descriptor array for current substance
        report(spin_system,'    building basis set descriptor...');
        basis_in_subst=cell(size(subgraphs_in_subst,1),1);
        parfor n=1:size(subgraphs_in_subst,1)

            % Isolate current subgraph
            subgraph=subgraphs_in_subst(n,:);
        
            % Determine the number of states in the
            % complete basis set of the current subgraph
            nstates_in_subgraph=prod(spin_dims_in_subst(subgraph)); %#ok<PFBNS>
        
            % Index and count the spins in the current subgraph
            spins_in_subgraph=find(subgraph); nspins_in_subgraph=nnz(subgraph);

            % Preallocate subgraph basis descriptor array 
            subgraph_basis=zeros(nstates_in_subgraph,...
                                 nspins_in_subgraph,state_idx_type);
        
            % Populate the local descriptor array
            for k=1:nspins_in_subgraph
            
                % Compute preceding dimension
                prefix_dim=prod(spin_dims_in_subst(spins_in_subgraph(1:(k-1))));
            
                % Get the current spin states
                current_spin_states=spin_state_lists_in_subst{spins_in_subgraph(k)}; %#ok<PFBNS>
            
                % Compute subsequent dimension
                suffix_dim=prod(spin_dims_in_subst(spins_in_subgraph((k+1):end)));
            
                % Combinatorial merge of preceding, current, and subsequent spin state lists
                subgraph_basis(:,k)=repmat(repelem(current_spin_states,suffix_dim,1),prefix_dim,1);

            end
        
            % Apply coherence order filter
            if ~isempty(projections_subst)

                % Compute coherence orders by an explicit loop to save memory
                coherence_orders=zeros([nstates_in_subgraph 1],state_idx_type);
                for k=1:nstates_in_subgraph

                    % Get projection quantum numbers
                    [~,M]=lin2lm(subgraph_basis(k,:));

                    % Compute coherence order
                    coherence_orders(k)=sum(M);

                end

                % Unit state is always to be kept
                keep_mask=false(nstates_in_subgraph,1); 
                keep_mask(1)=true();
            
                % Keep specified states
                for k=projections_subst
                    keep_mask=keep_mask|(coherence_orders==k);
                end
            
                % Remove the undesired states
                subgraph_basis(~keep_mask,:)=[];
                nstates_in_subgraph=size(subgraph_basis,1);

            end

            % Apply state space restriction 
            % by zero-quantum subspaces
            if any(zq_spins_in_subst,'all')
                
                % See if we have pertinent spins here
                pertinent_spins=zq_spins_in_subst(subgraph);
                if any(pertinent_spins,'all')

                    % Get subset coherence orders by an explicit loop to save memory
                    subset_coherence_orders=zeros([nstates_in_subgraph 1],state_idx_type);
                    for k=1:nstates_in_subgraph
                        [~,M]=lin2lm(subgraph_basis(k,pertinent_spins));
                        subset_coherence_orders(k)=sum(M);
                    end

                    % Unit state is always to be kept
                    keep_mask=false(nstates_in_subgraph,1);
                    keep_mask(1)=true();

                    % Keep ZQ states in the subset
                    keep_mask=keep_mask|(subset_coherence_orders==0);

                    % Remove the undesired states
                    subgraph_basis(~keep_mask,:)=[];
                    nstates_in_subgraph=size(subgraph_basis,1);

                end
               
            end

            % Return the indexing from subgraph into substance basis
            basis_in_subst{n}=spalloc(nstates_in_subgraph,nspins_in_subst,...
                                      numel(subgraph_basis),'single');
            basis_in_subst{n}(:,spins_in_subgraph)=subgraph_basis;

            % Release the memory
            subgraph_basis=[]; %#ok<NASGU>

        end
            
        % Pull basis descriptors from workers
        basis_in_subst=vertcat(basis_in_subst{:});
       
        % Fast repetition elimination using a hash table
        report(spin_system,'    eliminating redundant states...');
        basis_in_subst=unihash(basis_in_subst);
    
        % Lexicographic sorting of the basis
        report(spin_system,'    sorting the basis...');
        if (~isworkernode)&&(nnz(basis_in_subst)>1e5)

            % Run multi-threaded sorting
            basis_in_subst=distrib_dim(basis_in_subst,2);
            basis_in_subst=sortrows(basis_in_subst);
            basis_in_subst=gather(basis_in_subst);

        else

            % Run sorting in a single thread
            basis_in_subst=sortrows(basis_in_subst);

        end

        % Report the final number of states
        nstates_in_subst=size(basis_in_subst,1);
        report(spin_system,['    final number of basis states: ' ...
                                 int2str(nstates_in_subst)]);

        % Run rank and projection analysis 
        % in an explicit loop to save memory
        tot_proj=zeros([nstates_in_subst 1],state_idx_type);
        tot_cord=zeros([nstates_in_subst 1],state_idx_type);
        report(spin_system,'    computing quantum number indices...');
        parfor n=1:nstates_in_subst
            [L,M]=lin2lm(basis_in_subst(n,:));
            tot_proj(n)=sum(M,2);
            tot_cord(n)=sum(logical(L),2);
        end
        
        % Assign to the global object
        spin_system.bas.basis{s}=basis_in_subst; 
        spin_system.bas.tot_proj{s}=tot_proj;
        spin_system.bas.tot_cord{s}=tot_cord;

        % Clear huge variables
        clear('basis_in_subst','tot_proj','tot_cord');

    end

    % Process Zeeman basis sets
    if ismember(spin_system.bas.formalism,{'zeeman-wavef','zeeman-hilb','zeeman-liouv'})

        % Preallocate basis set array for the current substance,
        % same data type is safe because Hilbert space is smaller
        basis_in_subst=zeros(prod(mults_in_subst),...
                             nspins_in_subst,state_idx_type);

        % Fill basis set array
        for n=1:nspins_in_subst

            % Compute preceding dimension
            prefix_dim=prod(mults_in_subst(1:(n-1)));

            % Sequentially number the energy levels
            current_spin_states=1:mults_in_subst(n);
            current_spin_states=cast(current_spin_states,state_idx_type);
            current_spin_states=current_spin_states(:);

            % Compute subsequent dimension
            suffix_dim=prod(mults_in_subst((n+1):end));

            % Combinatorial merge of preceding, current, and subsequent energy level lists
            basis_in_subst(:,n)=repmat(repelem(current_spin_states,suffix_dim,1),prefix_dim,1);

        end

        % Store the dimension for diagnostics
        basis_dim_in_subst=size(basis_in_subst,1);

        % Report to the user
        switch spin_system.bas.formalism

            case 'zeeman-wavef'

                % Operators and wavefunctions
                report(spin_system,['    operator dimension:     '       ...
                                         int2str(basis_dim_in_subst) 'x' ...
                                         int2str(basis_dim_in_subst)]);
                report(spin_system,['    state vector dimension: '       ...
                                         int2str(basis_dim_in_subst) 'x1']);
                report(spin_system,' ');

            case 'zeeman-hilb'

                % Operators and density matrices
                report(spin_system,['    operator dimension:       '     ...
                                         int2str(basis_dim_in_subst) 'x' ...
                                         int2str(basis_dim_in_subst)]);
                report(spin_system,['    density matrix dimension: '     ...
                                         int2str(basis_dim_in_subst) 'x' ...
                                         int2str(basis_dim_in_subst)]);
                report(spin_system,' ');

            case 'zeeman-liouv'

                % Superoperators and state vectors
                report(spin_system,['    operator dimension:       '       ...
                                         int2str(basis_dim_in_subst^2) 'x' ...
                                         int2str(basis_dim_in_subst^2)]);
                report(spin_system,['    state vector dimension: '         ...
                                         int2str(basis_dim_in_subst^2) 'x1']);
                report(spin_system,' ');

        end
        
        % Store in the global object
        spin_system.bas.basis{s}=basis_in_subst;

    end

end

% Store the state count for each chemical substance
spin_system.bas.nstates=cellfun(@(x)size(x,1),spin_system.bas.basis);

% Spherical tensor basis has useful summaries 
if strcmp(spin_system.bas.formalism,'sphten-liouv')
    summary(spin_system,'basis_summary');
end

% Run the symmetry treatment
spin_system=symmetry(spin_system);

% Preload Lie algebra structure tables into RAM
if strcmp(spin_system.bas.formalism,'sphten-liouv')

    % Inform the user
    report(spin_system,'caching Lie structure tables...');

    % Find the spin multiplicities present
    unique_mults=unique(spin_system.comp.mults);

    % Preallocate the structure table arrays
    spin_system.bas.lpst=cell(max(unique_mults),1);
    spin_system.bas.rpst=cell(max(unique_mults),1);

    % Fill the arrays
    for n=setdiff(unique_mults,1)
        
        % Load from disk or compute
        [lpst,rpst]=ist_product_table(n);

        % Left product structure table
        spin_system.bas.lpst{n}=lpst;

        % Right product structure table
        spin_system.bas.rpst{n}=rpst;

    end

end

% Hash basis descriptor for caching tools later
if ismember('op_cache',spin_system.sys.enable)||...
   ismember('ham_cache',spin_system.sys.enable)
    report(spin_system,'computing basis set array hash...');
    spin_system.bas.basis_hash=md5_hash(spin_system.bas.basis);
end

end

% Consistency enforcement
function grumble(spin_system,bas)

% Check bas.formalism
if ~isfield(bas,'formalism')
    error('basis specification in bas.formalism is required.');
elseif ~ischar(bas.formalism)
    error('bas.formalism must be a string.');
elseif ~ismember(bas.formalism,{'zeeman-hilb','zeeman-liouv',...
                                'sphten-liouv','zeeman-wavef'})
    error('unrecognized formalism - see the basis preparation section of the manual.');
end

% Check zeeman-hilb formalism options
if strcmp(bas.formalism,'zeeman-hilb')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for zeeman-hilb formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'none'})
        error('bas.approximation should be set to ''none'' in zeeman-hilb formalism.');
    end

end

% Check zeeman-liouv formalism options
if strcmp(bas.formalism,'zeeman-liouv')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for zeeman-liouv formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'none'})
        error('bas.approximation should be set to ''none'' in zeeman-liouv formalism.');
    end

end

% Check sphten-liouv formalism options
if strcmp(bas.formalism,'sphten-liouv')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for sphten-liouv formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','IK-DNP','none'})
        error('unrecognized approximation - see the basis preparation section of the manual.');
    end
    
    % Check bas.connectivity
    if ismember(bas.approximation,{'IK-1','IK-2'})
        if ~isfield(bas,'connectivity')
            error('connectivity type must be specified in bas.connectivity variable.');
        elseif ~ischar(bas.connectivity)
            error('bas.connectivity must be a string.');
        elseif ~ismember(bas.connectivity,{'scalar_couplings','full_tensors'})
            error('unknown connectivity type - see the basis preparation section of the manual.');
        end
    end
    
    % Check bas.lever_level
    if ismember(bas.approximation,{'IK-0','IK-1','IK-DNP'})&&(~isfield(bas,'inter_level'))
        error('connectivity tracing depth must be specified in bas.inter_level variable.');
    end
    if ismember(bas.approximation,{'IK-0','IK-1'})
        if (~isnumeric(bas.inter_level))||(~isscalar(bas.inter_level))||...
           (mod(bas.inter_level,1)~=0)||(bas.inter_level<1)
            error('bas.inter_level must be a positive integer.');
        end
        if bas.inter_level>numel(spin_system.comp.isotopes)
            error('bas.inter_level cannot be greater than the number of spins in the system.');
        end
    end
    if strcmp(bas.approximation,'IK-DNP')
        if (~isnumeric(bas.inter_level))||(numel(bas.inter_level)~=3)||...
           any(mod(bas.inter_level,1)~=0,'all')||any(bas.inter_level<1,'all')
            error('bas.inter_level must be a vector with three positive integers.');
        end
        n_electrons=nnz(cellfun(@iselectron,spin_system.comp.isotopes));
        n_nuclei=nnz(cellfun(@isnucleus,spin_system.comp.isotopes));
        n_spins=numel(spin_system.comp.isotopes);
        if bas.inter_level(1)>n_electrons
            error('bas.inter_level(1) cannot exceed the number of electrons in the system.');
        end
        if bas.inter_level(2)>n_spins
            error('bas.inter_level(2) cannot exceed the number of spins in the system.');
        end
        if bas.inter_level(3)>n_nuclei
            error('bas.inter_level(3) cannot exceed the number of nuclei in the system.');
        end
    end
    
    % Check bas.prox_level
    if ismember(bas.approximation,{'IK-1','IK-2'})&&(~isfield(bas,'prox_level'))
        error('proximity tracing depth must be specified in bas.prox_level variable.');
    end
    if isfield(bas,'space_level')
        if  (~isnumeric(bas.prox_level))||(~isscalar(bas.prox_level))||...
            (mod(bas.prox_level,1)~=0)||(bas.prox_level<1)
            error('bas.prox_level must be a positive integer.');
        end
        if bas.prox_level>numel(spin_system.comp.isotopes)
            error('bas.prox_level cannot be greater than the number of spins in the system.');
        end
    end
    
    % Check bas.manual
    if isfield(bas,'manual')
        if (~islogical(bas.manual))&&(~isnumeric(bas.manual))
            error('bas.manual must be a logical matrix.');
        elseif size(bas.manual,2)~=spin_system.comp.nspins
            error('the number of rows in bas.manual must be equal to the number of spins in the system.')
        end
    end
    
    % Check bas.projections
    if isfield(bas,'projections')
        if ~iscell(bas.projections)
            error('bas.projections must be a cell array of vectors.');
        end
        for n=1:numel(bas.projections)
            if (~isnumeric(bas.projections{n}))||...
               (~isrow(bas.projections{n}))||any(mod(bas.projections{n},1)~=0)
                error('elements of bas.projections must be row vectors of integers.');
            end
        end
    end
    
    % Check bas.longitudinal
    if isfield(bas,'longitudinal')
        if ~iscell(bas.longitudinal)
            error('bas.longitudinal must be a cell array.');
        end
        for n=1:numel(bas.longitudinal)
            if ~iscell(bas.longitudinal{n})
                error('elements of bas.longitudinal must be a cell arrays.');
            end
            for k=1:numel(bas.longitudinal{n})
                if isnumeric(bas.longitudinal{n}{k})
                    if (~isreal(bas.longitudinal{n}{k}))||...
                       (~isscalar(bas.longitudinal{n}{k}))||...
                       (bas.longitudinal{n}{k}<1)||(mod(bas.longitudinal{n}{k},1)~=0)
                        error('numerical elements of bas.longitudinal must be positive integers.');
                    end
                    if bas.longitudinal{n}{k}>spin_system.comp.nspins
                        error('an element of bas.longitudinal exceeds the number of spins in the system.');
                    end
                elseif ischar(bas.longitudinal{n}{k})
                    if ~ismember(bas.longitudinal{n}{k},spin_system.comp.isotopes)
                        error('an element of bas.longitudinal refers to spins that are not present.');
                    end
                else
                    error('unrecognised bas.longitudinal specification.');
                end
            end
        end
    end

    % Check bas.zero_quantum
    if isfield(bas,'zero_quantum')
        if ~iscell(bas.zero_quantum)
            error('bas.zero_quantum must be a cell array.');
        end
        for n=1:numel(bas.zero_quantum)
            if ~iscell(bas.zero_quantum{n})
                error('elements of bas.zero_quantum must be a cell arrays.');
            end
            for k=1:numel(bas.zero_quantum{n})
                if isnumeric(bas.zero_quantum{n}{k})
                    if (~isreal(bas.zero_quantum{n}{k}))||...
                       (~isscalar(bas.zero_quantum{n}{k}))||...
                       (bas.zero_quantum{n}{k}<1)||(mod(bas.zero_quantum{n}{k},1)~=0)
                        error('numerical elements of bas.zero_quantum must be positive integers.');
                    end
                    if bas.zero_quantum{n}{k}>spin_system.comp.nspins
                        error('an element of bas.zero_quantum exceeds the number of spins in the system.');
                    end
                elseif ischar(bas.zero_quantum{n}{k})
                    if ~ismember(bas.zero_quantum{n}{k},spin_system.comp.isotopes)
                        error('an element of bas.zero_quantum refers to spins that are not present.');
                    end
                else
                    error('unrecognised bas.zero_quantum specification.');
                end
            end
        end
    end
    
end
    
% Disallow inapplicable approximations
if isfield(bas,'level')
    if ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','IK-DNP'})
        error('bas.inter_level is only applicable to IK-0,1,2,DNP basis sets.');
    end
end
if isfield(bas,'space_level')
    if ~ismember(bas.approximation,{'IK-1','IK-2'})
        error('bas.prox_level is only applicable to IK-1,2 basis sets.');
    end
end
if isfield(bas,'connectivity')
    if ~ismember(bas.approximation,{'IK-1','IK-2'})
        error('bas.connectivity is only applicable to IK-1,2 basis sets.');
    end
end

% Enforce sphten-liouv with criterion-based state pre-selection
if isfield(bas,'projections')&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('bas.projections option is only available for sphten-liouv formalism.');
end
if isfield(bas,'longitudinal')&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('bas.longitudinal option is only available for sphten-liouv formalism.');
end
if isfield(bas,'zero_quantum')&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('bas.zero_quantum option is only available for sphten-liouv formalism.');
end

% Enforce sphten-liouv when any kind of chemistry is present
if (numel(spin_system.chem.parts)>1)||(~isempty(spin_system.chem.flux_rate))
    if ~strcmp(bas.formalism,'sphten-liouv')
        error('chemical reaction modelling is only available for sphten-liouv formalism.');
    end
end

end

% In 1969, Robert Rathbun Wilson, the US physicist who headed Fermilab, the world's
% highest-energy particle accelerator laboratory, addressed the Congressional Joint
% Committee on Atomic Energy. Rhode Island Senator John Pastore asked Wilson to spell
% out what research into high-energy particle physics would do to improve the defence
% of the United States. Wilson gave a reply that went down in scientific history. Fer-
% milab, he said, had "nothing to do directly with defending our country, except to
% make it worth defending".
%
% http://www.theregister.co.uk/2009/02/09/woudhuysen_energise_1/

