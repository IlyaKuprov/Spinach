% Basis set control. This is the second mandatory function (after create.m)
% that must be called to build spin_system data structure. Syntax:
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
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=basis.m>

function spin_system=basis(spin_system,bas)

% Show the banner
banner(spin_system,'basis_banner');

% Check the input
grumble(spin_system,bas);

% Store the settings
spin_system.bas=bas;

% Find electrons and nuclei
e_idx=cellfun(@iselectron,spin_system.comp.isotopes);
n_idx=cellfun(@isnucleus,spin_system.comp.isotopes);

% Find bosonic modes
b_idx=ismember(spin_system.comp.types,{'C','V','T'});

% Report back to the user
summary_basis_opts(spin_system);

% Remind the user about the amplitude cut-off
report(spin_system,['coupling tensors with norm below ' ...
                    num2str(spin_system.tols.inter_cutoff) ...
                    ' Hz will be ignored.']);

% Process spherical tensor basis sets
if strcmp(spin_system.bas.formalism,'sphten-liouv')

    % Disallow spherical tensor basis sets for large multiplicities
    if any(spin_system.comp.mults>16,'all')
        error('multiplicities above 16 are not supported by sphten-liouv formalism.');
    end

    % Count chemical substances
    nsubst=numel(spin_system.chem.parts);

    % Run connectivity analysis for IK-DNP basis set
    if strcmp(spin_system.bas.approximation,'IK-DNP')

        % Make sure there are only electrons and nuclei
        if (nnz(e_idx)==0)||(nnz(n_idx)==0)
            error('IK-DNP approximation requires both electrons and nuclei.');
        end
        if ~all(e_idx|n_idx)
            error('IK-DNP approximation can only handle electrons and nuclei.');
        end

        % Isolate inter-electron interactions
        ee_couplings=spin_system.inter.coupling.matrix;
        ee_couplings(:,n_idx)={[]}; 
        ee_couplings(n_idx,:)={[]};

        % Isolate electron-nuclear interactions
        en_couplings=spin_system.inter.coupling.matrix;
        en_couplings(e_idx,e_idx)={[]}; 
        en_couplings(n_idx,n_idx)={[]};
        
        % Isolate inter-nuclear interactions
        nn_couplings=spin_system.inter.coupling.matrix;
        nn_couplings(:,e_idx)={[]}; 
        nn_couplings(e_idx,:)={[]};         

        % Remind the user about the amplitude cut-off
        report(spin_system,['coupling tensors with norm below ' ...
                            num2str(spin_system.tols.inter_cutoff) ...
                            ' Hz will be ignored.']);

        % Generate three types of connectivity graphs (e-e, e-n, n-n)
        ee_conmatrix=sparse(cellfun(@(x)norm(x,2),ee_couplings)>2*pi*spin_system.tols.inter_cutoff);
        en_conmatrix=sparse(cellfun(@(x)norm(x,2),en_couplings)>2*pi*spin_system.tols.inter_cutoff);
        nn_conmatrix=sparse(cellfun(@(x)norm(x,2),nn_couplings)>2*pi*spin_system.tols.inter_cutoff);

        % Make sure each spin is connected to itself
        ee_conmatrix=ee_conmatrix|speye(size(ee_conmatrix));
        en_conmatrix=en_conmatrix|speye(size(en_conmatrix));
        nn_conmatrix=nn_conmatrix|speye(size(nn_conmatrix));

        % Make sure connectivity is reciprocal
        ee_conmatrix=ee_conmatrix|transpose(ee_conmatrix);
        en_conmatrix=en_conmatrix|transpose(en_conmatrix);
        nn_conmatrix=nn_conmatrix|transpose(nn_conmatrix);

    end

    % Coupling tensor norm for the connectivity analysis
    if ismember(spin_system.bas.approximation,{'IK-1','IK-2','IK-SBS'})
        switch bas.connectivity

            case 'scalar_couplings'

                % Update the user
                report(spin_system,'scalar couplings will be used to build the coupling graph.');

                % Norm of the isotropic part of an interaction tensor
                tensor_norm=@(x)abs(trace(x)/3);

            case 'full_tensors'

                % Update the user
                report(spin_system,'full coupling tensors will be used to build the coupling graph.');

                % Full norm of an interaction tensor
                tensor_norm=@(x)norm(x,2);

        end
    end

    % Run connectivity analysis for IK-1,2 basis sets
    if ismember(spin_system.bas.approximation,{'IK-1','IK-2'})

        % Build the spin-spin coupling graph connectivity matrix
        inter_norm=cellfun(tensor_norm,spin_system.inter.coupling.matrix);
        spin_system.inter.conmatrix=sparse(inter_norm>2*pi*spin_system.tols.inter_cutoff);

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

    % Run connectivity analysis for IK-SBS basis set
    if strcmp(spin_system.bas.approximation,'IK-SBS')

        % Make sure there are both spins and bosonic modes
        if (nnz(b_idx)==0)||(nnz((~b_idx)&(spin_system.comp.mults>1))==0)
            error('IK-SBS approximation requires both spins and bosonic modes.');
        end

        % Build the spin-spin coupling graph connectivity matrix
        inter_norm=cellfun(tensor_norm,spin_system.inter.coupling.matrix);
        spin_system.inter.conmatrix=sparse(inter_norm>2*pi*spin_system.tols.inter_cutoff);

        % Bosonic mode connectivity
        if isfield(spin_system.inter,'modes')

            % Pairwise bosonic interaction fields list
            chan_flds={'exchange','dispersive','kerr','longitudinal'};

            % Largest pairwise coupling norm between modes and between modes and spins
            mode_norm=zeros(spin_system.comp.nspins);
            for n=1:numel(chan_flds)
                mode_norm=max(mode_norm,cellfun(@(x)norm(x,2),spin_system.inter.modes.(chan_flds{n})));
            end

            % Build the mode connectivity matrix
            mode_conmat=mode_norm>2*pi*spin_system.tols.inter_cutoff;

            % Over coupling tensor derivatives with respect to mode coordinates
            [m1,m2]=find(~cellfun(@isempty,spin_system.inter.modes.coupling_mod));
            for n=1:numel(m1)

                % Over derivative orders
                derivs=spin_system.inter.modes.coupling_mod{m1(n),m2(n)};
                for m=1:numel(derivs)

                    % Skip empty derivative orders
                    if isempty(derivs{m}), continue; end

                    % Spin pairs whose coupling derivative is above the cut-off
                    [p,q]=find(cellfun(tensor_norm,derivs{m})>2*pi*spin_system.tols.inter_cutoff);

                    % Connect the modes and the spin pair
                    for k=1:numel(p)
                        particles=[m1(n) m2(n) p(k) q(k)];
                        mode_conmat(particles,particles)=true;
                    end

                end

            end

            % Over effective field derivatives with respect to mode coordinates
            [m1,m2]=find(~cellfun(@isempty,spin_system.inter.modes.zeeman_mod));
            for n=1:numel(m1)

                % Over derivative orders
                derivs=spin_system.inter.modes.zeeman_mod{m1(n),m2(n)};
                for m=1:numel(derivs)

                    % Skip empty derivative orders
                    if isempty(derivs{m}), continue; end

                    % Spins whose field derivative is above the cut-off
                    p=find(cellfun(@(x)norm(x,2),derivs{m})>2*pi*spin_system.tols.inter_cutoff);

                    % Connect the modes and the spin
                    for k=1:numel(p)
                        particles=[m1(n) m2(n) p(k)];
                        mode_conmat(particles,particles)=true;
                    end

                end

            end

            % Merge spin and bosonic mode connectivity matrices
            spin_system.inter.conmatrix=spin_system.inter.conmatrix|sparse(mode_conmat);

            % Update the user
            report(spin_system,'bosonic mode couplings above the cut-off added to the coupling graph.');

        end

        % Make sure connectivity is reciprocal
        spin_system.inter.conmatrix=spin_system.inter.conmatrix|transpose(spin_system.inter.conmatrix);

        % Isolate boson-boson connectivity
        bb_conmatrix=spin_system.inter.conmatrix;
        bb_conmatrix(~b_idx,:)=false;
        bb_conmatrix(:,~b_idx)=false;

        % Isolate spin-boson connectivity
        sb_conmatrix=spin_system.inter.conmatrix;
        sb_conmatrix(b_idx,b_idx)=false;
        sb_conmatrix(~b_idx,~b_idx)=false;

        % Isolate spin-spin connectivity
        ss_conmatrix=spin_system.inter.conmatrix;
        ss_conmatrix(b_idx,:)=false;
        ss_conmatrix(:,b_idx)=false;

        % Make sure each particle is connected to itself
        bb_conmatrix=bb_conmatrix|speye(size(bb_conmatrix));
        sb_conmatrix=sb_conmatrix|speye(size(sb_conmatrix));
        ss_conmatrix=ss_conmatrix|speye(size(ss_conmatrix));

        % Issue a report to the user
        report(spin_system,['boson-boson connectivity matrix density ' num2str(100*nnz(bb_conmatrix)/numel(bb_conmatrix)) '%']);
        report(spin_system,['spin-boson connectivity matrix density ' num2str(100*nnz(sb_conmatrix)/numel(sb_conmatrix)) '%']);
        report(spin_system,['spin-spin connectivity matrix density ' num2str(100*nnz(ss_conmatrix)/numel(ss_conmatrix)) '%']);

    end

    % Smallest signed integer class that holds every single-spin state index
    idx_class=min_int_type(max([1 spin_system.comp.mults.^2-1]),'signed');

    % Build state lists for individual spins
    spin_state_lists=cell(spin_system.comp.nspins,1);
    for n=1:spin_system.comp.nspins
        spin_state_lists{n}=cast((0:(spin_system.comp.mults(n)^2-1))',idx_class);
    end

    % Apply longitudinal filters
    if isfield(bas,'longitudinal')
        for s=1:nsubst
            for k=1:numel(bas.longitudinal{s})

                % Find the specified spins within the current substance
                if isnumeric(bas.longitudinal{s}{k})
                    spins_in_question=bas.longitudinal{s}{k}(:)';
                    if ~all(ismember(spins_in_question,spin_system.chem.parts{s}))
                        error(['bas.longitudinal{' int2str(s) '} refers to spins outside substance ' int2str(s) '.']);
                    end
                else
                    spins_in_question=spin_system.chem.parts{s}(strcmp(bas.longitudinal{s}{k},...
                                      spin_system.comp.isotopes(spin_system.chem.parts{s})));
                    if isempty(spins_in_question)
                        error(['no ' bas.longitudinal{s}{k} ' spins in substance ' int2str(s) '.']);
                    end
                end

                % Kill unwanted states
                for n=spins_in_question(:)'
                    report(spin_system,['keeping only longitudinal states on spin ' num2str(n) '...']);
                    [~,M]=lin2lm(spin_state_lists{n}); spin_state_lists{n}(M~=0)=[];
                end

            end
        end
    end

    % Compute subspace dimensions for individual spins
    spin_dims=cellfun(@numel,spin_state_lists);

    % Preallocate subgraph lists and their substance indices
    subgraphs=cell(nsubst,1); subgraph_subst=cell(nsubst,1);

    % Loop over chemical substances
    for s=1:nsubst

        % Spins of the current substance
        spins_in_subst=spin_system.chem.parts{s}(:)';
        nspins_in_subst=numel(spins_in_subst);
        report(spin_system,['chemical substance ' int2str(s) ', ' int2str(nspins_in_subst) ' spins:']);

        % Generate subgraphs in the spin index of the substance
        switch spin_system.bas.approximation

            case 'none'

                % Single subgraph with all spins of the substance
                coupling_subgraphs=true(1,nspins_in_subst);

                % Do not run proximity analysis
                proximity_subgraphs=false(0,nspins_in_subst);

            case 'IK-0'

                % Clip the correlation level to the substance size
                inter_level=min([nspins_in_subst bas.inter_level]);

                % Find all possible groups of inter_level spins
                col_index=nchoosek(1:nspins_in_subst,inter_level);

                % Get the number of groups
                ngroups=size(col_index,1);

                % Assign numbers to groups
                row_index=repmat((1:ngroups)',1,inter_level);

                % Generate the subgraph list
                coupling_subgraphs=sparse(row_index,col_index,true,ngroups,nspins_in_subst);
                report(spin_system,['    ' num2str(ngroups) ' subgraphs generated by combinatorial analysis.']);

                % Do not run proximity analysis
                proximity_subgraphs=false(0,nspins_in_subst);

            case 'IK-1'

                % Clip the correlation levels to the substance size
                inter_level=min([nspins_in_subst bas.inter_level]);
                prox_level=min([nspins_in_subst bas.prox_level]);

                % Run connectivity analysis
                coupling_subgraphs=dfpt(spin_system.inter.conmatrix(spins_in_subst,spins_in_subst),inter_level);
                report(spin_system,['    ' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);

                % Run proximity analysis
                proximity_subgraphs=dfpt(spin_system.inter.proxmatrix(spins_in_subst,spins_in_subst),prox_level);
                report(spin_system,['    ' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);

            case 'IK-2'

                % Clip the proximity level to the substance size
                prox_level=min([nspins_in_subst bas.prox_level]);

                % Run connectivity analysis
                coupling_subgraphs=unique(spin_system.inter.conmatrix(spins_in_subst,spins_in_subst),'rows');
                report(spin_system,['    ' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);

                % Run proximity analysis
                proximity_subgraphs=dfpt(spin_system.inter.proxmatrix(spins_in_subst,spins_in_subst),prox_level);
                report(spin_system,['    ' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);

            case 'IK-DNP'

                % Clip the correlation levels to the substance size
                inter_level=min(nspins_in_subst,bas.inter_level);

                % Inter-electron connectivity analysis
                ee_subgraphs=dfpt(ee_conmatrix(spins_in_subst,spins_in_subst),inter_level(1));
                report(spin_system,['    generated ' num2str(size(ee_subgraphs,1)-nnz(n_idx(spins_in_subst))) ' inter-electron subgraphs.']);

                % Electron-nuclear connectivity analysis
                en_subgraphs=dfpt(en_conmatrix(spins_in_subst,spins_in_subst),inter_level(2));
                report(spin_system,['    generated ' num2str(size(en_subgraphs,1)) ' electron-nuclear subgraphs.']);

                % Inter-nuclear connectivity analysis
                nn_subgraphs=dfpt(nn_conmatrix(spins_in_subst,spins_in_subst),inter_level(3));
                report(spin_system,['    generated ' num2str(size(nn_subgraphs,1)-nnz(e_idx(spins_in_subst))) ' inter-nuclear subgraphs.']);

                % Merge coupling subgraph lists
                coupling_subgraphs=[ee_subgraphs; en_subgraphs; nn_subgraphs];

                % Do not run proximity analysis
                proximity_subgraphs=false(0,nspins_in_subst);

            case 'IK-SBS'

                % Clip the correlation levels to the substance size
                inter_level=min(nspins_in_subst,bas.inter_level);

                % Boson-boson connectivity analysis
                bb_subgraphs=dfpt(bb_conmatrix(spins_in_subst,spins_in_subst),inter_level(1));
                report(spin_system,['    generated ' num2str(size(bb_subgraphs,1)-nnz(~b_idx(spins_in_subst))) ' boson-boson subgraphs.']);

                % Spin-boson connectivity analysis
                sb_subgraphs=dfpt(sb_conmatrix(spins_in_subst,spins_in_subst),inter_level(2));
                report(spin_system,['    generated ' num2str(size(sb_subgraphs,1)) ' spin-boson subgraphs.']);

                % Spin-spin connectivity analysis
                ss_subgraphs=dfpt(ss_conmatrix(spins_in_subst,spins_in_subst),inter_level(3));
                report(spin_system,['    generated ' num2str(size(ss_subgraphs,1)-nnz(b_idx(spins_in_subst))) ' spin-spin subgraphs.']);

                % Merge coupling subgraph lists
                coupling_subgraphs=[bb_subgraphs; sb_subgraphs; ss_subgraphs];

                % Do not run proximity analysis
                proximity_subgraphs=false(0,nspins_in_subst);

            otherwise

                % Complain and bomb out
                error('unrecognised basis set.');

        end

        % Include user-specified subgraphs that belong to the substance
        manual_subgraphs=false(0,nspins_in_subst);
        if isfield(bas,'manual')
            for n=1:size(bas.manual,1)
                if nnz(bas.manual(n,spins_in_subst))==nnz(bas.manual(n,:))
                    manual_subgraphs=[manual_subgraphs; logical(bas.manual(n,spins_in_subst))]; %#ok<AGROW>
                elseif nnz(bas.manual(n,spins_in_subst))>0
                    error(['row ' int2str(n) ' of bas.manual crosses chemical substance boundaries.']);
                end
            end
            report(spin_system,['    added ' num2str(size(manual_subgraphs,1)) ' subgraphs specified by the user.']);
        end

        % Assemble the subgraph list
        subgraphs_in_subst=[coupling_subgraphs; proximity_subgraphs; manual_subgraphs];
        clear('coupling_subgraphs','proximity_subgraphs','manual_subgraphs');

        % Prune subgraphs involving spin zero particles
        subgraphs_in_subst(:,spin_system.comp.mults(spins_in_subst)==1)=false;

        % Remove empty, identical, and enclosed subgraphs
        subgraphs_in_subst=subgraphs_in_subst(any(subgraphs_in_subst,2),:);
        subgraphs_in_subst=prune_subgraphs(unique(subgraphs_in_subst,'rows'));

        % Report back to the user
        subgraph_sizes=sum(subgraphs_in_subst,2);
        for n=min(subgraph_sizes):max(subgraph_sizes)
            if nnz(subgraph_sizes==n)>0
                report(spin_system,['    keeping ' num2str(nnz(subgraph_sizes==n)) ' subgraphs with ' num2str(n) ' spins each.']);
            end
        end
        if isfield(bas,'projections')&&(~isempty(bas.projections{s}))
            report(spin_system,['    keeping only coherence orders with M=[' num2str(bas.projections{s}) ']...']);
        end

        % Embed the subgraphs into the full spin index
        subgraphs{s}=false(size(subgraphs_in_subst,1),spin_system.comp.nspins);
        subgraphs{s}(:,spins_in_subst)=subgraphs_in_subst;
        subgraph_subst{s}=repmat(s,[size(subgraphs_in_subst,1) 1]);

    end

    % Merge the subgraph lists of all substances
    subgraphs=vertcat(subgraphs{:}); subgraph_subst=vertcat(subgraph_subst{:});
    clear('subgraphs_in_subst','subgraph_sizes');

    % Resolve zero-quantum filter spins for each substance
    zq_spins=false(nsubst,spin_system.comp.nspins);
    if isfield(bas,'zero_quantum')
        for s=1:nsubst
            for k=1:numel(bas.zero_quantum{s})

                % Find the specified spins within the current substance
                if isnumeric(bas.zero_quantum{s}{k})
                    spins_in_question=bas.zero_quantum{s}{k}(:)';
                    if ~all(ismember(spins_in_question,spin_system.chem.parts{s}))
                        error(['bas.zero_quantum{' int2str(s) '} refers to spins outside substance ' int2str(s) '.']);
                    end
                else
                    spins_in_question=spin_system.chem.parts{s}(strcmp(bas.zero_quantum{s}{k},...
                                      spin_system.comp.isotopes(spin_system.chem.parts{s})));
                    if isempty(spins_in_question)
                        error(['no ' bas.zero_quantum{s}{k} ' spins in substance ' int2str(s) '.']);
                    end
                end

                % Add to the zero-quantum spin set of the substance
                zq_spins(s,spins_in_question)=true();

            end
            if any(zq_spins(s,:))
                report(spin_system,['keeping only the zero-quantum states on spins ' num2str(find(zq_spins(s,:))) '...']);
            end
        end
    end

    % Balance the subgraph list
    shuffle=randperm(size(subgraphs,1));
    subgraphs=subgraphs(shuffle,:); subgraph_subst=subgraph_subst(shuffle);

    % Populate the basis descriptor array
    report(spin_system,'building basis set descriptor...');
    basis_spec=cell(size(subgraphs,1),1);
    parfor n=1:size(subgraphs,1)

        % Determine which spins belong to the current subgraph
        spins_involved=find(subgraphs(n,:)); nspins_involved=numel(spins_involved);

        % Determine the total number of states in the current subgraph
        nstates=prod(spin_dims(spins_involved)); %#ok<PFBNS>

        % Preallocate the local descriptor array
        local_basis_spec=zeros(nstates,nspins_involved,idx_class);

        % Populate the local descriptor array
        for k=1:nspins_involved

            % Compute preceding dimension
            dim_before=prod(spin_dims(spins_involved(1:(k-1))));

            % Get the current spin states
            current_states=spin_state_lists{spins_involved(k)}; %#ok<PFBNS>

            % Compute following dimension
            dim_after=prod(spin_dims(spins_involved((k+1):end)));

            % Replicate the state list into the direct product order
            local_basis_spec(:,k)=repmat(repelem(current_states,dim_after,1),dim_before,1);

        end

        % Apply coherence order filter, always keeping the unit state
        if isfield(bas,'projections')&&(~isempty(bas.projections{subgraph_subst(n)}))
            [~,M]=lin2lm(local_basis_spec);
            state_mask=ismember(sum(M,2),bas.projections{subgraph_subst(n)});
            state_mask(1)=true(); local_basis_spec=local_basis_spec(state_mask,:);
        end

        % Apply zero-quantum filter over the specified spins of the substance
        if any(zq_spins(subgraph_subst(n),spins_involved)) %#ok<PFBNS>
            [~,M]=lin2lm(local_basis_spec(:,zq_spins(subgraph_subst(n),spins_involved)));
            local_basis_spec=local_basis_spec(sum(M,2)==0,:);
        end

        % Drop excessive inter-nuclear correlations
        if strcmp(spin_system.bas.approximation,'IK-DNP') %#ok<PFBNS>

            % Identify inter-nuclear correlations
            corr_order=sum(logical(local_basis_spec),2);
            nn_corr_order=sum(logical(local_basis_spec(:,n_idx(spins_involved))),2); %#ok<PFBNS>
            pure_nn_state=(corr_order==nn_corr_order);
            pure_nn_state(1)=false;

            % Build the drop mask
            state_mask=pure_nn_state&(nn_corr_order>bas.inter_level(3));

            % Kill the states
            local_basis_spec(state_mask,:)=[];

        end

        % Drop excessive pure boson-boson and pure spin-spin correlations
        if strcmp(spin_system.bas.approximation,'IK-SBS')

            % Identify pure boson-boson and pure spin-spin correlations
            corr_order=sum(logical(local_basis_spec),2);
            bb_corr_order=sum(logical(local_basis_spec(:,b_idx(spins_involved))),2); %#ok<PFBNS>
            ss_corr_order=corr_order-bb_corr_order;
            pure_bb_state=(corr_order==bb_corr_order);
            pure_ss_state=(corr_order==ss_corr_order);
            pure_bb_state(1)=false; pure_ss_state(1)=false;

            % Build the drop mask
            state_mask=(pure_bb_state&(bb_corr_order>bas.inter_level(1)))|...
                       (pure_ss_state&(ss_corr_order>bas.inter_level(3)));

            % Kill the states
            local_basis_spec(state_mask,:)=[];

        end

        % Embed the descriptor into the full spin index
        [rows,cols,vals]=find(local_basis_spec);
        basis_spec{n}=sparse(rows,spins_involved(cols),double(vals),size(local_basis_spec,1),spin_system.comp.nspins);

    end

    % Deallocate variables
    clear('spin_state_lists','subgraphs','subgraph_subst','spin_dims','zq_spins');

    % Pull basis descriptor from the nodes, unit state first
    basis_spec=[sparse(1,spin_system.comp.nspins); vertcat(basis_spec{:})];

    % Eliminate redundant states and sort, on a dense integer copy if that is smaller
    report(spin_system,'eliminating redundant states and sorting the basis...');
    idx_bytes=numel(typecast(zeros(1,idx_class),'uint8'));
    if numel(basis_spec)*idx_bytes<=16*nnz(basis_spec)+8*(spin_system.comp.nspins+1)
        [rows,cols,vals]=find(basis_spec);
        basis_spec=zeros(size(basis_spec),idx_class);
        basis_spec(sub2ind(size(basis_spec),rows,cols))=vals;
        basis_spec=unique(basis_spec,'rows');
        [rows,cols,vals]=find(basis_spec);
        spin_system.bas.basis=sparse(rows,cols,double(vals),size(basis_spec,1),spin_system.comp.nspins);
    else
        spin_system.bas.basis=unique(basis_spec,'rows');
    end

    % Deallocate variables
    clear('basis_spec');

    % Total projection quantum number and correlation order of each state
    [L,M]=lin2lm(spin_system.bas.basis);
    spin_system.bas.tot_proj=full(sum(M,2));
    spin_system.bas.tot_cord=full(sum(logical(L),2));
    clear('L','M');

    % Report on chemical species
    for s=1:nsubst
        nstates=nnz(any(spin_system.bas.basis(:,spin_system.chem.parts{s}),2));
        report(spin_system,['chemical substance ' int2str(s) ': ' num2str(nstates) ' states.']);
    end

    % Print the summary
    summary_basis(spin_system);

    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);

end

% Process Hilbert space Zeeman basis
if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-wavef'})

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

% Process Liouville space Zeeman basis
if strcmp(spin_system.bas.formalism,'zeeman-liouv')

    % Build the Hilbert space Zeeman index table
    dim=prod(spin_system.comp.mults);
    zbas=zeros(dim,spin_system.comp.nspins);
    for n=1:spin_system.comp.nspins
        current_column=1;
        for k=1:spin_system.comp.nspins
            if n==k
                current_column=kron(current_column,(1:spin_system.comp.mults(k))');
            else
                current_column=kron(current_column,ones(spin_system.comp.mults(k),1));
            end
        end
        zbas(:,n)=current_column;
    end

    % Ket and bra index tables in the vectorisation order
    spin_system.bas.basis=[repmat(zbas,[dim 1]) kron(zbas,ones(dim,1))];

    % Report to the user
    report(spin_system,['matrix dimension for all superoperators and state vectors: ' num2str(dim^2)]);

    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);

end

% Preload Lie algebra structure tables
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

% Hash the basis descriptor for caching tools later
if ismember('op_cache',spin_system.sys.enable)||...
   ismember('ham_cache',spin_system.sys.enable)
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
    elseif ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','IK-DNP','IK-SBS','none'})
        error('unrecognized approximation - see the basis preparation section of the manual.');
    end

    % Disallow bosonic modes in IK-1,2 basis sets
    if ismember(bas.approximation,{'IK-1','IK-2'})&&any(ismember(spin_system.comp.types,{'C','V','T'}))
        error('IK-1 and IK-2 basis sets are for spin-only systems, use IK-SBS when bosonic modes are present.');
    end

    % Check bas.connectivity
    if ismember(bas.approximation,{'IK-1','IK-2','IK-SBS'})
        if ~isfield(bas,'connectivity')
            error('connectivity type must be specified in bas.connectivity variable.');
        elseif ~ischar(bas.connectivity)
            error('bas.connectivity must be a string.');
        elseif ~ismember(bas.connectivity,{'scalar_couplings','full_tensors'})
            error('unknown connectivity type - see the basis preparation section of the manual.');
        end
    end

    % Check bas.inter_level
    if ismember(bas.approximation,{'IK-0','IK-1','IK-DNP','IK-SBS'})&&(~isfield(bas,'inter_level'))
        error('connectivity tracing depth must be specified in bas.inter_level variable.');
    end
    if ismember(bas.approximation,{'IK-0','IK-1'})
        if (~isnumeric(bas.inter_level))||(~isscalar(bas.inter_level))||(mod(bas.inter_level,1)~=0)||(bas.inter_level<1)
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
    if strcmp(bas.approximation,'IK-SBS')
        if (~isnumeric(bas.inter_level))||(numel(bas.inter_level)~=3)||...
           any(mod(bas.inter_level,1)~=0,'all')||any(bas.inter_level<1,'all')
            error('bas.inter_level must be a vector with three positive integers.');
        end
        mode_mask=ismember(spin_system.comp.types,{'C','V','T'});
        n_modes=nnz(mode_mask); n_spins=nnz((~mode_mask)&(spin_system.comp.mults>1));
        if bas.inter_level(1)>n_modes
            error('bas.inter_level(1) cannot exceed the number of bosonic modes in the system.');
        end
        if bas.inter_level(2)>n_modes+n_spins
            error('bas.inter_level(2) cannot exceed the number of particles in the system.');
        end
        if bas.inter_level(3)>n_spins
            error('bas.inter_level(3) cannot exceed the number of spins in the system.');
        end
    end

    % Check bas.prox_level
    if ismember(bas.approximation,{'IK-1','IK-2'})&&(~isfield(bas,'prox_level'))
        error('proximity tracing depth must be specified in bas.prox_level variable.');
    end
    if isfield(bas,'prox_level')
        if  (~isnumeric(bas.prox_level))||(~isscalar(bas.prox_level))||(mod(bas.prox_level,1)~=0)||(bas.prox_level<1)
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
            error('the number of columns in bas.manual must be equal to the number of spins in the system.');
        end
    end

    % Check bas.projections
    if isfield(bas,'projections')
        if (~iscell(bas.projections))||(numel(bas.projections)~=numel(spin_system.chem.parts))
            error('bas.projections must be a cell array with one element per chemical substance.');
        end
        for n=1:numel(bas.projections)
            if (~isnumeric(bas.projections{n}))||(~isreal(bas.projections{n}))||...
               ((~isempty(bas.projections{n}))&&(~isrow(bas.projections{n})))||...
               any(mod(bas.projections{n},1)~=0,'all')
                error('elements of bas.projections must be empty or row vectors of integers.');
            end
        end
    end

    % Check bas.longitudinal
    if isfield(bas,'longitudinal')
        if (~iscell(bas.longitudinal))||(numel(bas.longitudinal)~=numel(spin_system.chem.parts))
            error('bas.longitudinal must be a cell array with one element per chemical substance.');
        end
        for n=1:numel(bas.longitudinal)
            if ~iscell(bas.longitudinal{n})
                error('elements of bas.longitudinal must be cell arrays.');
            end
            for k=1:numel(bas.longitudinal{n})
                if isnumeric(bas.longitudinal{n}{k})
                    if (~isreal(bas.longitudinal{n}{k}))||...
                       any(mod(bas.longitudinal{n}{k},1)~=0,'all')||...
                       any(bas.longitudinal{n}{k}<1,'all')||...
                       any(bas.longitudinal{n}{k}>spin_system.comp.nspins,'all')
                        error('numeric entries in bas.longitudinal must be positive integers within the system bounds.');
                    end
                elseif ischar(bas.longitudinal{n}{k})
                    if ~ismember(bas.longitudinal{n}{k},spin_system.comp.isotopes)
                        error('bas.longitudinal refers to spins that are not present in the system.');
                    end
                else
                    error('bas.longitudinal must contain isotope strings or vectors of spin numbers.');
                end
            end
        end
    end

    % Check bas.zero_quantum
    if isfield(bas,'zero_quantum')
        if (~iscell(bas.zero_quantum))||(numel(bas.zero_quantum)~=numel(spin_system.chem.parts))
            error('bas.zero_quantum must be a cell array with one element per chemical substance.');
        end
        for n=1:numel(bas.zero_quantum)
            if ~iscell(bas.zero_quantum{n})
                error('elements of bas.zero_quantum must be cell arrays.');
            end
            for k=1:numel(bas.zero_quantum{n})
                if isnumeric(bas.zero_quantum{n}{k})
                    if (~isreal(bas.zero_quantum{n}{k}))||...
                       any(mod(bas.zero_quantum{n}{k},1)~=0,'all')||...
                       any(bas.zero_quantum{n}{k}<1,'all')||...
                       any(bas.zero_quantum{n}{k}>spin_system.comp.nspins,'all')
                        error('numeric entries in bas.zero_quantum must be positive integers within the system bounds.');
                    end
                elseif ischar(bas.zero_quantum{n}{k})
                    if ~ismember(bas.zero_quantum{n}{k},spin_system.comp.isotopes)
                        error('bas.zero_quantum refers to spins that are not present in the system.');
                    end
                else
                    error('bas.zero_quantum must contain isotope strings or vectors of spin numbers.');
                end
            end
        end
    end

end

% Catch retired field names
if isfield(bas,'level')
    error('bas.level has been renamed bas.inter_level.');
end
if isfield(bas,'space_level')
    error('bas.space_level has been renamed bas.prox_level.');
end
if isfield(bas,'longitudinals')
    error('bas.longitudinals has been renamed bas.longitudinal, with one cell array per chemical substance.');
end

% Disallow inapplicable approximations
if isfield(bas,'inter_level')
    if ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','IK-DNP','IK-SBS'})
        error('bas.inter_level is only applicable to IK-0,1,2,DNP,SBS basis sets.');
    end
end
if isfield(bas,'prox_level')
    if ~ismember(bas.approximation,{'IK-1','IK-2'})
        error('bas.prox_level is only applicable to IK-1,2 basis sets.');
    end
end
if isfield(bas,'connectivity')
    if ~ismember(bas.approximation,{'IK-1','IK-2','IK-SBS'})
        error('bas.connectivity is only applicable to IK-1,2,SBS basis sets.');
    end
end

% Enforce sphten-liouv with state filters
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

