% Extended validation of user-declared permutation symmetry. Confirms
% that the Zeeman, coupling, and giant-spin interactions stored in the
% spin system object are strictly invariant under every operation of
% each declared permutation group, so that the irreducible representa-
% tion projectors built by symmetry.m correspond to a symmetry that
% the interactions actually possess. The declared symmetry is a permu-
% tation of spin labels, not a spatial rotation; interaction tensors
% related by a rotation rather than being identical are not accepted.
% When no symmetry is declared, the function returns without perfor-
% ming any checks. Syntax:
%
%                    validate_sym(spin_system,bas)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object, with the
%                   interaction arrays already processed by create.m
%
%    bas          - basis specification structure; the fields used
%                   here are bas.sym_group (cell array of group names)
%                   and bas.sym_spins (cell array of spin index vec-
%                   tors), as described in the basis specification
%                   section of the online manual
%
% Outputs:
%
%    this function returns nothing; it throws a descriptive error when
%    the interaction data does not obey a declared symmetry
%
% Note: this is a service function of the Spinach kernel that should
%       not be called directly; it is called by symmetry.m
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=validate_sym.m>

function validate_sym(spin_system,bas)

% Check consistency
grumble(spin_system,bas);

% Skip validation when symmetry is switched off
if ismember('symmetry',spin_system.sys.disable), return; end

% Skip validation when no symmetry group is declared
if (~isfield(bas,'sym_group'))||isempty(bas.sym_group), return; end

% Interaction agreement tolerance, rad/s
tol=2*pi*spin_system.tols.inter_cutoff;

% Number of spins in the system
nspins=spin_system.comp.nspins;

% Zeeman and giant-spin interaction arrays
zeeman=spin_system.inter.zeeman.matrix;
giant=spin_system.inter.giant.coeff;

% Loop over declared symmetry groups
for m=1:numel(bas.sym_group)

    % Spins in the current symmetry group
    spins=bas.sym_spins{m};

    % All spins in a symmetry group must be the same isotope
    if any(~strcmp(spin_system.comp.isotopes(spins),spin_system.comp.isotopes{spins(1)}))
        error(['spins declared symmetric under ' bas.sym_group{m} ...
               ' are not all of the same isotope.']);
    end

    % Permutation elements of the declared group
    group=perm_group(bas.sym_group{m});

    % Loop over symmetry operations
    for n=1:group.order

        % Global spin permutation for this operation
        perm=1:nspins; perm(spins)=spins(group.elements(n,:));

        % Loop over spins in the symmetry group
        for k=1:numel(spins)

            % Check Zeeman tensor invariance
            if norm(zeeman{perm(spins(k))}-zeeman{spins(k)},2)>tol
                error(['spins ' num2str(spins(k)) ' and ' num2str(perm(spins(k))) ...
                       ' are declared symmetric under ' bas.sym_group{m} ', but their '...
                       'Zeeman interaction tensors are not identical.']);
            end

            % Check giant-spin coefficient invariance across all ranks
            coeff_a=giant{spins(k)}; coeff_b=giant{perm(spins(k))};
            giant_mismatch=(numel(coeff_a)~=numel(coeff_b));
            for r=1:numel(coeff_a)
                if (~giant_mismatch)&&(norm(coeff_b{r}-coeff_a{r},2)>tol)
                    giant_mismatch=true;
                end
            end
            if giant_mismatch
                error(['spins ' num2str(spins(k)) ' and ' num2str(perm(spins(k))) ...
                       ' are declared symmetric under ' bas.sym_group{m} ', but their '...
                       'giant-spin interactions are not identical.']);
            end

            % Check coupling tensor invariance against every spin
            for q=1:nspins
                if norm(get_coupling(spin_system,perm(spins(k)),perm(q))-...
                        get_coupling(spin_system,spins(k),q),2)>tol % #NORMOK
                    error(['the interaction between spins ' num2str(spins(k)) ' and ' ...
                           num2str(q) ' is not preserved by the ' bas.sym_group{m} ...
                           ' symmetry; it differs from the interaction between spins ' ...
                           num2str(perm(spins(k))) ' and ' num2str(perm(q)) '; the condition '...
                           'enforced is componentwise identity in the laboratory frame, '...
                           'because symmetry.m permutes spin labels without rotating it.']);
                end
            end

        end

    end

end

% Confirm success to the user
report(spin_system,'declared permutation symmetry is consistent with the interaction data.');

end

% Consistency enforcement
function grumble(spin_system,bas)

if isfield(bas,'sym_group')
    if (~iscell(bas.sym_group))||any(~cellfun(@ischar,bas.sym_group))
        error('bas.sym_group must be a cell array of strings.');
    end
    if ~isfield(bas,'sym_spins')
        error('bas.sym_spins must be specified alongside bas.sym_group.');
    end
    if (~iscell(bas.sym_spins))||any(~cellfun(@isnumeric,bas.sym_spins))
        error('bas.sym_spins must be a cell array of vectors.');
    end
    if numel(bas.sym_spins)~=numel(bas.sym_group)
        error('bas.sym_group and bas.sym_spins arrays must have the same number of elements.');
    end
    for n=1:numel(bas.sym_spins)
        if any(bas.sym_spins{n}>spin_system.comp.nspins)||any(bas.sym_spins{n}<1)||(numel(bas.sym_spins{n})<2)
            error('incorrect spin labels in bas.sym_spins.');
        end
    end
end

end

