% Basis set control. This is the second mandatory function (after create.m)
% that must be called in every calculation to build spin_system data struc-
% ture. Syntax:
%
%                     spin_system=basis(spin_system,bas)
%
% Parameters:
%
%     spin_system   - primary Spinach data structure, the output
%                     of create.m function
%
%     bas           - basis set specification structure described
%                     in detail in the online manual
%
% Outputs:
%
%     spin_system   - primary Spinach data structure, updated with
%                     the basis set and related information
%
% Note: it is important to understand the factors that influence basis set
%       selection in spin dynamics simulations - see our paper
%
%                http://link.aip.org/link/doi/10.1063/1.3624564
%
%       for further information on this subject.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=basis.m>

function spin_system=basis(spin_system,bas)

% Show the banner
banner(spin_system,'basis_banner');

% Check the input
grumble(spin_system,bas);

% Store the settings
spin_system.bas=bas;

% Report back to the user
summary(spin_system,'basis_settings');

% Process spherical tensor basis sets
if strcmp(spin_system.bas.formalism,'sphten-liouv')

    % Run connectivity analysis
    if ismember(spin_system.bas.approximation,{'IK-1','IK-2'})
        
        % Build connectivity and proximity matrices
        switch bas.connectivity
            
            case 'scalar_couplings'
                
                % Use scalar parts of all interaction tensors
                report(spin_system,'scalar couplings will be used to build the coupling graph.');
                spin_system.inter.conmatrix=sparse(abs(cellfun(@trace,spin_system.inter.coupling.matrix)/3)>2*pi*spin_system.tols.inter_cutoff);
                
            case 'full_tensors'
                
                % Use complete interaction tensors
                report(spin_system,'full coupling tensors will be used to build the coupling graph.');
                spin_system.inter.conmatrix=sparse(cellfun(@(x)norm(x,2),spin_system.inter.coupling.matrix)>2*pi*spin_system.tols.inter_cutoff);
                
        end
        
        % Remind the user about the amplitude cut-off
        report(spin_system,['coupling tensors with norm below ' num2str(spin_system.tols.inter_cutoff) ' Hz will be ignored.']);
        
        % Make sure each spin is connected and proximate to itself
        spin_system.inter.conmatrix=spin_system.inter.conmatrix|speye(size(spin_system.inter.conmatrix));
        spin_system.inter.proxmatrix=spin_system.inter.proxmatrix|speye(size(spin_system.inter.proxmatrix));
        
        % Make sure connectivity and proximity are reciprocal
        spin_system.inter.conmatrix=spin_system.inter.conmatrix|transpose(spin_system.inter.conmatrix);
        spin_system.inter.proxmatrix=spin_system.inter.proxmatrix|transpose(spin_system.inter.proxmatrix);
        
        % Issue a report to the user
        report(spin_system,['connectivity matrix density ' num2str(100*nnz(spin_system.inter.conmatrix)/numel(spin_system.inter.conmatrix)) '%']);
        report(spin_system,['proximity matrix density ' num2str(100*nnz(spin_system.inter.proxmatrix)/numel(spin_system.inter.proxmatrix)) '%']);
        
        % Determine the number of independent subsystems
        n_subsystems=max(scomponents(spin_system.inter.conmatrix|spin_system.inter.proxmatrix));
        
        % Print a notice to the user
        if n_subsystems>1
            report(spin_system,['WARNING - there are ' num2str(n_subsystems) ' subsystems that are not coupled to each other.']);
        end
        
    end
    
    % Build state lists for individual spins
    spin_state_lists=cell(spin_system.comp.nspins,1);
    for n=1:spin_system.comp.nspins
        spin_state_lists{n}=(0:(spin_system.comp.mults(n)^2-1))';
    end
    
    % Apply longitudinal filters
    if isfield(bas,'longitudinals')
        
        % Build the spin list
        spins_in_question=[];
        for n=1:numel(bas.longitudinals)
            if isnumeric(bas.longitudinals{n})
                spins_in_question=[spins_in_question bas.longitudinals{n}(:)']; %#ok<AGROW>
            else
                spins_in_question=[spins_in_question find(strcmp(bas.longitudinals{n},spin_system.comp.isotopes))]; %#ok<AGROW>
            end
        end
        spins_in_question=unique(spins_in_question(:))';
        
        % Kill unwanted states
        for n=spins_in_question
            report(spin_system,['keeping only longitudinal states on spin ' num2str(n) '...']);
            [~,M]=lin2lm(spin_state_lists{n}); spin_state_lists{n}(M~=0)=[];
        end
        
    end
    
    % Compute subspace dimensions for individual spins
    spin_dims=cellfun(@numel,spin_state_lists);

    % Generate subgraphs
    switch spin_system.bas.approximation 
        
        case 'none'
            
            % Match the chemical subsystems
            coupling_subgraphs=false(numel(spin_system.chem.parts),spin_system.comp.nspins);
            for n=1:numel(spin_system.chem.parts)
                coupling_subgraphs(n,spin_system.chem.parts{n})=true();
            end
            
            % Do not run proximity analysis
            proximity_subgraphs=[];
            
        case 'IK-0'
            
            % Make the coupling subgraphs array
            coupling_subgraphs=cell(numel(spin_system.chem.parts),1);
            
            % Loop over chemical subsystems
            for n=1:numel(spin_system.chem.parts)
                
                % Find all possible groups of bas.level spins
                col_index=nchoosek(spin_system.chem.parts{n},bas.level);
                
                % Get the number of groups
                ngroups=size(col_index,1);
                
                % Assign numbers to groups
                row_index=repmat((1:ngroups)',1,bas.level);
                
                % Generate the subgraph list
                coupling_subgraphs{n}=sparse(row_index,col_index,true,...
                                             ngroups,spin_system.comp.nspins);
                
            end
            
            % Merge coupling subgraph lists
            coupling_subgraphs=cell2mat(coupling_subgraphs);
            
            % Inform the user
            report(spin_system,[num2str(size(coupling_subgraphs,1)) ...
                                ' subgraphs generated by combinatorial analysis.']);
            
            % Do not run proximity analysis
            proximity_subgraphs=[];
        
        case 'IK-1'
        
            % Run connectivity analysis
            coupling_subgraphs=dfpt(spin_system.inter.conmatrix,bas.level);
            report(spin_system,['' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);
            
            % Run proximity analysis
            proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
            report(spin_system,['' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);
        
        case 'IK-2'
        
            % Run connectivity analysis
            coupling_subgraphs=spin_system.inter.conmatrix;
            report(spin_system,['' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);
            
            % Run proximity analysis
            proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
            report(spin_system,['' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);
            
        otherwise
            
            % Complain and bomb out
            error('unrecognised basis set.');
        
    end

    % Include user-specified subgraphs
    if (~isfield(bas,'manual'))||isempty(bas.manual)
        manual_subgraphs=[];
    else
        manual_subgraphs=bas.manual;
        report(spin_system,['added ' num2str(size(manual_subgraphs,1)) ' subgraphs specified by the user.']);
    end
    
    % Assemble the subgraph list
    subgraphs=[coupling_subgraphs; proximity_subgraphs; manual_subgraphs];
    clear('coupling_subgraphs','proximity_subgraphs','manual_subgraphs');
    
    % Prune subgraphs involving spin zero particles
    report(spin_system,'pruning subgraphs involving zero spin particles...');
    subgraphs(:,spin_system.comp.mults==1)=0;
    
    % Remove identical subgraphs
    report(spin_system,'removing identical subgraphs...');
    subgraphs=unique(subgraphs,'rows');
    
    % Store subgraphs for future use if needed
    if ismember('xmemlist',spin_system.sys.enable)
        spin_system.bas.subgraphs=logical(subgraphs);
    end
    
    % Report back to the user
    subgraph_sizes=sum(subgraphs,2);
    for n=min(subgraph_sizes):max(subgraph_sizes)
        if nnz(subgraph_sizes==n)>0
            report(spin_system,['generated ' num2str(nnz(subgraph_sizes==n)) ' subgraphs with ' num2str(n) ' spins each.']);
        end
    end
    if isfield(bas,'projections')
        report(spin_system,['keeping only coherence orders with M=[' num2str(bas.projections) ']...']);
    end
    if isfield(bas,'zero_quantum')
        for n=1:numel(bas.zero_quantum)
            if isnumeric(bas.zero_quantum{n})
                report(spin_system,['keeping only the zero-quantum states on spins ' num2str(bas.zero_quantum{n}) '...']);
            else
                report(spin_system,['keeping only the zero-quantum states on ' bas.zero_quantum{n} '...']);
            end
        end
    end
        
    % Balance the subgraph list
    subgraphs=subgraphs(randperm(size(subgraphs,1)),:);
    
    % Populate the basis descriptor array
    report(spin_system,'building basis set descriptor...');
    basis_spec=cell(size(subgraphs,1),1);
    parfor n=1:size(subgraphs,1)
        
        % Determine the total number of states in the current subgraph
        nstates=prod(spin_dims(logical(subgraphs(n,:)))); %#ok<PFBNS>
        
        % Determine which spins belong to the current subgraph
        spins_involved=find(subgraphs(n,:));
        
        % Preallocate the local descriptor array
        local_basis_spec=spalloc(nstates,spin_system.comp.nspins,nstates*nnz(subgraphs(n,:))); %#ok<PFBNS>
        
        % Populate the local descriptor array
        for k=1:numel(spins_involved)
            
            % Compute preceding dimension
            dim_before=prod(spin_dims(spins_involved(1:(k-1))));
            
            % Get the current spin states
            current_states=spin_state_lists{spins_involved(k)}; %#ok<PFBNS>
            
            % Compute following dimension
            dim_after=prod(spin_dims(spins_involved((k+1):end)));
            
            % Kron everything together
            local_basis_spec(:,spins_involved(k))=kron(kron(ones(dim_before,1),current_states),ones(dim_after,1)); %#ok<SPRIX>
            
        end
        
        % Apply coherence order filter
        if isfield(bas,'projections')
            
            % Compute coherence order for each basis element
            [~,M]=lin2lm(local_basis_spec);
            
            % Start with empty mask
            state_mask=false(size(local_basis_spec,1),1);
            
            % Keep specified coherence orders and the unit state
            state_mask(1)=true(); projection_numbers=sum(M,2);
            for k=bas.projections
                state_mask=state_mask|(projection_numbers==k);
            end
            
            % Kill the undesired states
            local_basis_spec(~state_mask,:)=[]; %#ok<SPRIX>
            
        end
        
        % Apply zero-quantum filter
        if isfield(bas,'zero_quantum')
            
            % Start with empty mask
            state_mask=false(size(local_basis_spec,1),1);
            
            % Loop over the specified spins
            for k=1:numel(bas.zero_quantum)
                
                % Find the specified spins
                if isnumeric(bas.zero_quantum{k})
                    spins_in_question=bas.zero_quantum{k};
                else
                    spins_in_question=strcmp(bas.zero_quantum{k},spin_system.comp.isotopes);
                end
                
                % Analyze the basis
                [~,M]=lin2lm(local_basis_spec(:,spins_in_question));
                
                % Update the mask
                state_mask=or(state_mask,sum(M,2)~=0);
                
            end
            
            % Kill the states
            local_basis_spec(state_mask,:)=[]; %#ok<SPRIX>
            
        end
        
        % Assign the global variable
        basis_spec{n}=local_basis_spec;
        
    end
    
    % Deallocate variables
    clear('spin_state_lists','local_basis_spec','subgraphs','spin_dims',...
          'current_states','subgraph_sizes','local_basis_hash');
    if isfield(bas,'projections')
        clear('state_mask','M','projection_numbers');
    end
    
    % Pull basis descriptor from the nodes
    basis_spec=vertcat(basis_spec{:});

    % Build a hash table
    basis_hash=repmat(' ',[size(basis_spec,1) 32]);
    parfor k=1:size(basis_spec,1)
        basis_hash(k,:)=md5_hash(full(basis_spec(k,:)));
    end
   
    % Eliminate redundant states using hash table
    report(spin_system,'eliminating redundant states...');
    [~,idx]=unique(basis_hash,'rows','stable'); 
    basis_spec=basis_spec(idx,:);
    
    % Deallocate variables
    clear('basis_hash', 'idx');
    
    % Sort the basis explicitly
    report(spin_system,'sorting the basis...');
    if (~isworkernode)&&(nnz(basis_spec)>1e5)
        
        % Run multithreaded sorting
        basis_spec=distrib_dim(basis_spec,2);
        basis_spec=sortrows(basis_spec);
        spin_system.bas.basis=gather(basis_spec);
        
    else
        
        % Run sorting in a single thread
        spin_system.bas.basis=sortrows(basis_spec);
        
    end
    
    % Deallocate variables
    clear('basis_spec');
   
    % Report on chemical species
    chem_idx=false(size(spin_system.bas.basis,1),numel(spin_system.chem.parts));
    for n=1:numel(spin_system.chem.parts)
        chem_idx(:,n)=(sum(spin_system.bas.basis(:,spin_system.chem.parts{n}),2)>0);
        report(spin_system,['chemical species ' num2str(n) ': ' num2str(nnz(chem_idx(:,n))) ' states.']);
    end
    
    % Make sure chemical species are unlinked
    if any(sum(chem_idx(2:end,:),2)~=1)
        error('some basis set elements belong to either none or multiple chemical species.');
    end
    
    % Build state-cluster cross-membership list if needed
    if ismember('xmemlist',spin_system.sys.enable)
        
        % Report to the user
        report(spin_system,'building state-subgraph cross-membership list... ');
        
        % Localise variables
        basis_loc=spin_system.bas.basis;
        subgraphs_loc=spin_system.bas.subgraphs;
        
        % Preallocate the result
        nstates=size(basis_loc,1); nclusters=size(subgraphs_loc,1);
        xmemlist=spalloc(nstates,nclusters,ceil(nstates*nclusters/spin_system.comp.nspins));
        
        % Run the matching
        parfor n=1:size(subgraphs_loc,1)
            xmemlist(:,n)=~any(basis_loc(:,~subgraphs_loc(n,:)),2); %#ok<SPRIX,PFBNS>
        end
        
        % Store the result
        spin_system.bas.xmemlist=xmemlist;

    end 
    
    % Print the summary
    summary(spin_system,'basis');

    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);

end

% Process Hilbert space Zeeman basis
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
   
    % Preallocate basis set array
    spin_system.bas.basis=zeros(prod(spin_system.comp.mults),spin_system.comp.nspins);
    
    % Fill basis set array
    for n=1:spin_system.comp.nspins
        current_column=1;
        for k=1:spin_system.comp.nspins
            if n==k
                current_column=kron(current_column,(1:spin_system.comp.mults(k))');
            else
                current_column=kron(current_column,ones(spin_system.comp.mults(k),1));
            end
        end
        spin_system.bas.basis(:,n)=current_column;
    end
    
    % Report to the user
    report(spin_system,['matrix dimension for all operators and states: ' num2str(prod(spin_system.comp.mults))]);
    
    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);
    
end

% Hash the basis descriptor for caching tools later
if ismember('op_cache',spin_system.sys.enable)
    spin_system.bas.basis_hash=md5_hash(spin_system.bas.basis);
end

end

% Grumble function
function grumble(spin_system,bas)

% Check bas.formalism
if ~isfield(bas,'formalism')
    error('basis specification in bas.formalism is required.');
elseif ~ischar(bas.formalism)
    error('bas.formalism must be a string.');
elseif ~ismember(bas.formalism,{'zeeman-hilb','zeeman-liouv','sphten-liouv'})
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
    elseif ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','none'})
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
    
    % Check bas.level
    if ismember(bas.approximation,{'IK-0','IK-1'})&&(~isfield(bas,'level'))
        error('connectivity tracing depth must be specified in bas.level variable.');
    end
    if isfield(bas,'level')
        if (~isnumeric(bas.level))||(~isscalar(bas.level))||(mod(bas.level,1)~=0)||(bas.level<1)
            error('bas.level must be a positive integer.');
        end
        if bas.level>numel(spin_system.comp.isotopes)
            error('bas.level cannot be greater than the number of spins in the system.');
        end
    end
    
    % Check bas.space_level
    if ismember(bas.approximation,{'IK-1','IK-2'})&&(~isfield(bas,'space_level'))
        error('proximity tracing depth must be specified in bas.space_level variable.');
    end
    if isfield(bas,'space_level')
        if  (~isnumeric(bas.space_level))||(~isscalar(bas.space_level))||(mod(bas.space_level,1)~=0)||(bas.space_level<1)
            error('bas.space_level must be a positive integer.');
        end
        if bas.space_level>numel(spin_system.comp.isotopes)
            error('bas.space_level cannot be greater than the number of spins in the system.');
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
        if (~isnumeric(bas.projections))||(~isrow(bas.projections))||any(mod(bas.projections,1)~=0)
            error('bas.projections must be a row vector of integers.');
        end
    end
    
    % Check bas.longitudinals
    if isfield(bas,'longitudinals')
        if (~iscell(bas.longitudinals))||any(~cellfun(@ischar,bas.longitudinals))
            error('bas.longitudinals must be a cell array of strings.');
        end
        if any(~ismember(bas.longitudinals,spin_system.comp.isotopes))
            error('bas.longitudinals refers to spins that are not present in the system.');
        end
    end
    
end
    
% Disallow inapplicable approximations
if isfield(bas,'level')
    if ~ismember(bas.approximation,{'IK-0','IK-1','IK-2'})
        error('bas.level is only applicable to IK-0,1,2 basis sets.');
    end
end
if isfield(bas,'space_level')
    if ~ismember(bas.approximation,{'IK-1','IK-2'})
        error('bas.space_level is only applicable to IK-1,2 basis sets.');
    end
end
if isfield(bas,'connectivity')
    if ~ismember(bas.approximation,{'IK-1','IK-2'})
        error('bas.connectivity is only applicable to IK-1,2 basis sets.');
    end
end

% Enforce sphten-liouv with projection selection
if isfield(bas,'projections')&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('bas.projections option is only available for sphten-liouv formalism.');
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

