% Extended validation of user-declared permutation symmetry. Confirms
% that the Zeeman and coupling interaction tensors stored in the spin
% system object are invariant under every operation of each declared
% symmetry group, so that the irreducible representation projectors
% built by symmetry.m correspond to a symmetry that the interactions
% actually possess. When no symmetry is declared, the function returns
% without performing any checks. Syntax:
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

% Zeeman and coupling interaction arrays
zeeman=spin_system.inter.zeeman.matrix;
coupling=spin_system.inter.coupling.matrix;

% Loop over declared symmetry groups
for m=1:numel(bas.sym_group)

    % Spins in the current symmetry group
    spins=bas.sym_spins{m};

    % Permutation elements of the declared group
    group=perm_group(bas.sym_group{m});

    % Loop over symmetry operations
    for n=1:group.order

        % Global spin permutation for this operation
        perm=1:nspins; perm(spins)=spins(group.elements(n,:));

        % Loop over spins in the symmetry group
        for k=1:numel(spins)

            % Check Zeeman tensor invariance under the permutation
            if norm(zeeman{perm(spins(k))}-zeeman{spins(k)},2)>tol
                error(['spins ' num2str(spins(k)) ' and ' num2str(perm(spins(k))) ...
                       ' are declared symmetric under ' bas.sym_group{m} ', but their '...
                       'Zeeman interaction tensors are not identical.']);
            end

            % Check interaction tensor invariance against every spin
            for q=1:nspins
                if norm(pair_coupling(coupling,perm(spins(k)),perm(q))-...
                        pair_coupling(coupling,spins(k),q),2)>tol
                    error(['the interaction between spins ' num2str(spins(k)) ' and ' ...
                           num2str(q) ' is not preserved by the ' bas.sym_group{m} ...
                           ' symmetry; it differs from the interaction between spins ' ...
                           num2str(perm(spins(k))) ' and ' num2str(perm(q)) '.']);
                end
            end

        end

    end

end

% Confirm success to the user
report(spin_system,'declared permutation symmetry is consistent with the interaction data.');

end

% Interaction tensor acting on spins p and q, with empty blocks read as zeros
function A=pair_coupling(coupling,p,q)

if ~isempty(coupling{p,q})
    A=coupling{p,q};
elseif ~isempty(coupling{q,p})
    A=coupling{q,p}';
else
    A=zeros(3);
end

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

