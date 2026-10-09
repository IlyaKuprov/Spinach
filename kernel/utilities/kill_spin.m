% Removes specified particles (spins or bosonic modes) from the
% spin_system structure and updates dependent data. Syntax:
%
%            spin_system=kill_spin(spin_system,hit_list)
%
% Parameters:
%
%     spin_system   - primary Spinach data structure
%
%     hit_list      - a vector of integers or a logical
%                     vector giving particle numbers in the
%                     unified isotope list to be removed
%
% Outputs:
%
%     spin_system   - the data structure with the indica-
%                     ted particles and dependent infor-
%                     mation updated; an existing basis is rebuilt
%
% Notes: an existing basis is rebuilt from its input settings, with
%        local manual columns and global filter labels reindexed.
%        Isotope filters with no surviving local spins are removed.
%        Retained depths are capped by the surviving local populations;
%        vector depths use the particle-type bounds enforced by basis.
%        IK-DNP or IK-SBS losing a required particle class switches to
%        IK-0 at the surviving class depth, without connectivity pruning.
%        Symmetry and assumption information is cleared; call assume
%        again before constructing a Hamiltonian. Mode strengths are
%        cleared; the mode container is removed when no bosonic
%        particles remain. Spin-free substances are retained.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=kill_spin.m>

function spin_system=kill_spin(spin_system,hit_list)

% Check consistency
grumble(spin_system,hit_list)

% Catch logical indexing
if islogical(hit_list), hit_list=find(hit_list); end

% Retain input settings rather than patching compiled basis data
if isfield(spin_system,'bas')
    fields={'formalism','approximation','inter_level','prox_level',...
            'space_level','connectivity','manual','projections',...
            'longitudinal','zero_quantum'};
    fields=intersect(fields,fieldnames(spin_system.bas));
    for n=1:numel(fields)
        bas.(fields{n})=spin_system.bas.(fields{n});
    end

    % Reindex retained filter labels and remove empty isotope selections
    keep=setdiff(1:spin_system.comp.nspins,hit_list);
    for n=1:numel(spin_system.chem.parts)
        local_keep=~ismember(spin_system.chem.parts{n},hit_list);
        if isfield(bas,'manual')
            bas.manual{n}=bas.manual{n}(:,local_keep);
        end
        for field={'longitudinal','zero_quantum'}
            if isfield(bas,field{1})
                for k=1:numel(bas.(field{1}){n})
                    labels=bas.(field{1}){n}{k};
                    if isnumeric(labels)
                        [present,labels]=ismember(labels,keep);
                        bas.(field{1}){n}{k}=labels(present);
                    elseif ~ismember(labels,spin_system.comp.isotopes(...
                                     spin_system.chem.parts{n}(local_keep)))
                        bas.(field{1}){n}{k}=[];
                    end
                end
                bas.(field{1}){n}=bas.(field{1}){n}(~cellfun(@isempty,bas.(field{1}){n}));
            end
        end

        % A substance losing its last spin retains only its unit coordinate
        if ~any(local_keep)
            bas.approximation{n}='none';
            for field={'inter_level','prox_level','space_level','connectivity'}
                if isfield(bas,field{1}), bas.(field{1}){n}=[]; end
            end
        else

            % Cap scalar depths by the surviving local particle count
            for field={'inter_level','prox_level','space_level'}
                if isfield(bas,field{1})
                    bas.(field{1}){n}=min(bas.(field{1}){n},nnz(local_keep));
                end
            end

            % Bound vector depths by the same populations used in basis validation
            spins=spin_system.chem.parts{n}(local_keep);
            if strcmp(bas.approximation{n},'IK-DNP')
                isotopes=spin_system.comp.isotopes(spins);
                bas.inter_level{n}(1)=min(bas.inter_level{n}(1),nnz(cellfun(@iselectron,isotopes)));
                bas.inter_level{n}(3)=min(bas.inter_level{n}(3),nnz(cellfun(@isnucleus,isotopes)));
            elseif strcmp(bas.approximation{n},'IK-SBS')
                modes=ismember(spin_system.comp.types(spins),{'C','V','T'});
                nspins=nnz(~modes&(spin_system.comp.mults(spins)>1));
                bas.inter_level{n}(1)=min(bas.inter_level{n}(1),nnz(modes));
                bas.inter_level{n}(2)=min(bas.inter_level{n}(2),nnz(modes)+nspins);
                bas.inter_level{n}(3)=min(bas.inter_level{n}(3),nspins);
            end

            % Retain the surviving class depth without a two-class graph requirement
            if ismember(bas.approximation{n},{'IK-DNP','IK-SBS'})&&...
               any(bas.inter_level{n}([1 3])==0)
                bas.approximation{n}='IK-0';
                bas.inter_level{n}=max(1,max(bas.inter_level{n}([1 3])));
                if isfield(bas,'connectivity'), bas.connectivity{n}=[]; end
            end
        end
    end
end

% Inform the user
report(spin_system,['removing ' num2str(numel(hit_list)) ...
                    ' particles from the system...']);

% Update isotope and particle type lists
spin_system.comp.isotopes(hit_list)=[];
spin_system.comp.types(hit_list)=[];

% Update the isotope list hash
if isfield(spin_system.comp,'iso_hash')
    spin_system.comp.iso_hash=md5_hash(spin_system.comp.isotopes);
end

% Update spin numbers
spin_system.comp.nspins=spin_system.comp.nspins-numel(hit_list);

% Update labels list
spin_system.comp.labels(hit_list)=[];

% Update multiplicities and magnetogyric ratios
spin_system.comp.mults(hit_list)=[];
spin_system.inter.gammas(hit_list)=[];

% Update base frequencies
spin_system.inter.basefrqs(hit_list)=[];

% Update Zeeman tensor array
spin_system.inter.zeeman.matrix(hit_list)=[];

% Update DD scaling multipliers
spin_system.inter.zeeman.ddscal(hit_list)=[];

% Update giant spin Hamiltonian terms
spin_system.inter.giant.coeff(hit_list)=[];

% Update coupling tensor array
spin_system.inter.coupling.matrix(hit_list,:)=[];
spin_system.inter.coupling.matrix(:,hit_list)=[];

% Update coordinates
spin_system.inter.coordinates(hit_list)=[];

% Update proximity matrix
spin_system.inter.proxmatrix(hit_list,:)=[];
spin_system.inter.proxmatrix(:,hit_list)=[];

% Update particle-indexed bosonic mode data
if isfield(spin_system.inter,'modes')

    % Remove particle coordinates from scalar mode parameters
    fields={'frqs','carriers','anharms','damp','dephase'};
    for n=1:numel(fields)
        spin_system.inter.modes.(fields{n})(hit_list)=[];
    end

    % Remove particle coordinates from mode pair channels
    fields={'exchange','kerr','longitudinal','dispersive',...
            'coupling_mod','zeeman_mod'};
    for n=1:numel(fields)
        pairs=spin_system.inter.modes.(fields{n});
        pairs(hit_list,:)=[]; pairs(:,hit_list)=[];

        % Reindex spin leaves inside retained modulation derivative orders
        if ismember(fields{n},{'coupling_mod','zeeman_mod'})
            for k=find(~cellfun(@isempty,pairs(:)))'
                orders=pairs{k};
                for p=1:numel(orders)
                    if isempty(orders{p}), continue; end
                    orders{p}(:,hit_list)=[];
                    if strcmp(fields{n},'coupling_mod')
                        orders{p}(hit_list,:)=[];
                    end
                end
                pairs{k}=orders;
            end
        end
        spin_system.inter.modes.(fields{n})=pairs;
    end

    % Discard inapplicable mode data and stale mode assumptions
    if ~any(ismember(spin_system.comp.types,{'C','V','T'}))
        spin_system.inter=rmfield(spin_system.inter,'modes');
    elseif isfield(spin_system.inter.modes,'strength')
        spin_system.inter.modes=rmfield(spin_system.inter.modes,'strength');
    end
end

% Update relaxation parameters
if ~isempty(spin_system.rlx.r1_rates)
    spin_system.rlx.r1_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.r2_rates)
    spin_system.rlx.r2_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.lind_r1_rates)
    spin_system.rlx.lind_r1_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.lind_r2_rates)
    spin_system.rlx.lind_r2_rates(hit_list)=[];
end
if ~isempty(spin_system.rlx.srfk_mdepth)
    spin_system.rlx.srfk_mdepth(hit_list,:)=[];
    spin_system.rlx.srfk_mdepth(:,hit_list)=[];
end
if ~isempty(spin_system.rlx.weiz_r1d)
    spin_system.rlx.weiz_r1d(hit_list,:)=[];
    spin_system.rlx.weiz_r1d(:,hit_list)=[];
end
if ~isempty(spin_system.rlx.weiz_r2d)
    spin_system.rlx.weiz_r2d(hit_list,:)=[];
    spin_system.rlx.weiz_r2d(:,hit_list)=[];
end

% Update scalar relaxation source spins
if ~isempty(spin_system.rlx.srsk_sources)
    srsk_spins=zeros(1,spin_system.comp.nspins+numel(hit_list));
    srsk_spins(spin_system.rlx.srsk_sources)=1; srsk_spins(hit_list)=[];
    spin_system.rlx.srsk_sources=find(srsk_spins);
end

% Update kinetics parameters
for n=1:numel(spin_system.chem.parts)
    subsystem_idx=false(1,spin_system.comp.nspins+numel(hit_list));
    subsystem_idx(spin_system.chem.parts{n})=true();
    subsystem_idx(hit_list)=[];
    spin_system.chem.parts{n}=find(subsystem_idx);
end

% Rebuild the reaction atom maps in the surviving global spin index
keep=setdiff(1:spin_system.comp.nspins+numel(hit_list),hit_list);
for n=1:numel(spin_system.chem.reactions)
    reaction=spin_system.chem.reactions{n};
    matching=reaction.matching;
    matching(any(ismember(matching,hit_list),2),:)=[];
    [~,reaction.matching]=ismember(matching,keep);
    if isfield(reaction,'selector')&&ischar(reaction.selector{1})
        [~,reaction.selector{2}]=ismember(reaction.selector{2},keep);
    end
    spin_system.chem.reactions{n}=reaction;
end

% If any basis set information is found, destroy it
if isfield(spin_system,'bas')
    spin_system=rmfield(spin_system,'bas');
end

% If any connectivity information is found, destroy it
if isfield(spin_system.inter,'conmatrix')
    spin_system.inter=rmfield(spin_system.inter,'conmatrix');
end

% If any symmetry information is found, destroy it
if isfield(spin_system.comp,'sym_group')
    spin_system.comp=rmfield(spin_system.comp,{'sym_group','sym_spins','sym_a1g_only'});
end

% If any assumption information is found, destroy it
if isfield(spin_system.inter,'assumptions')
    spin_system.inter=rmfield(spin_system.inter,'assumptions');
    report(spin_system,'WARNING - assumption information must be re-created.');
end
if isfield(spin_system.inter.zeeman,'strength')
    spin_system.inter.zeeman=rmfield(spin_system.inter.zeeman,'strength');
    report(spin_system,'WARNING - assumption information must be re-created.');
end
if isfield(spin_system.inter.giant,'strength')
    spin_system.inter.giant=rmfield(spin_system.inter.giant,'strength');
    report(spin_system,'WARNING - assumption information must be re-created.');
end
if isfield(spin_system.inter.coupling,'strength')
    spin_system.inter.coupling=rmfield(spin_system.inter.coupling,'strength');
    report(spin_system,'WARNING - assumption information must be re-created.');
end

% Rebuild descriptors, symmetry projectors, dimensions, offsets, and cache identity
if exist('bas','var'), spin_system=basis(spin_system,bas); end

end

% Consistency enforcement
function grumble(spin_system,hit_list)
if islogical(hit_list)
    if(numel(hit_list)~=spin_system.comp.nspins)
        error('the size of the hit mask does not match the number of spins in the system.');
    end
else
    if (~isnumeric(hit_list))||any(hit_list<1)||any(mod(hit_list,1)~=0)
        error('hit_list must be a logical mask or an array of positive integers.');
    end
    if any(hit_list>spin_system.comp.nspins)
        error('at least one number in hit_list exceeds the number of spins.');
    end
end
if islogical(hit_list), hit_list=find(hit_list); end
for n=1:numel(spin_system.chem.reactions)
    reaction=spin_system.chem.reactions{n};
    if isfield(reaction,'selector')
        if ischar(reaction.selector{1})
            if any(ismember(reaction.selector{2},hit_list))
                error('Spinach:kill_spin:selectorElectron','cannot remove an electron used by a reaction selector.');
            end
        elseif any(ismember(spin_system.chem.parts{reaction.reactants},hit_list))
            error('Spinach:kill_spin:selectorMatrix','rebuild user selector matrices before removing spins from their substance.');
        end
    end
end
end

% I do not see that the sex of the candidate is an argument against
% her admission as privatdozent. After all, we are a university, not
% a bath house.
%
% David Hilbert, about Emmy Noether, in 1915.

