% Finds out which substance hosts the specified spins;
% throws an error if there is more than one. Syntax:
%
%         subst=which_subst(spin_system,spins)
%
% Parameters:
%
%     spins - a list of positive integers
%
% Outputs:
%
%     subst - a positive integer
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=which_subst.m>

function subst=which_subst(spin_system,spins)

% Check consistency
grumble(spin_system,spins);

% Find the substances hosting specified spins
subst_mask=cellfun(@(x)any(ismember(spins,x)),...
                   spin_system.chem.parts);

% Only one substance is permitted
if nnz(subst_mask)>1
    error('spin list crosses chemical boundaries.');
elseif nnz(subst_mask)==0
    error('spins do not belong to any substance.');
end

% Get substance number
subst=find(subst_mask);

% Confirm that all spins are in the same substance
if ~all(ismember(spins,spin_system.chem.parts{subst}))
    error('spin list crosses chemical boundaries.');
end

end

% Consistency enforcement
function grumble(spin_system,spins)
if (~isnumeric(spins))||(~isreal(spins))||...
   (~isvector(spins))||any(mod(spins,1)~=0)||any(spins<1)
    error('spins must be an array of positive integers.');
end
if any(spins>spin_system.comp.nspins)
    error('spin number exceeds the total spin count.');
end
if numel(unique(spins))~=numel(spins)
    error('some spins are specified repeatedly.');
end
end

% I swear to you that to think too much is
% a disease. A real, actual disease.
%
% Fyodor Dostoevsky

