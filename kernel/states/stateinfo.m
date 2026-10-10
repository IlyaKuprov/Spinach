% Prints the state vector norm and the list of the most populated basis 
% states in the order of decreasing population. Syntax:
%
%                   stateinfo(spin_system,rho,npops)
%
% Parameters:
%
%                rho   - state vector
%
%              npops   - number of largest populations to print
%
% Outputs:
%
%   This function prints a summary of the state composition to the con-
%   sole in the following format:
%
%         (L1,M1) (L2,M2) ... (Ln,Mn)     coefficient     number
%
%   This corresponds to the direct product of single-spin irreducible
%   spherical tensors with the specified indices, its coefficient in 
%   the linear combination, and the number of the corresponding state
%   in the direct-sum basis set. Each row names its substance and
%   global spin indices; only that substance's tensor labels are shown.
%
% Note: this function requires a spherical tensor basis set.
%
% ilya.kuprov@weizmann.ac.il
% kpervushin@ntu.edu.sg
%
% <https://spindynamics.org/wiki/index.php?title=stateinfo.m>

function stateinfo(spin_system,rho,npops)

% Check consistency
grumble(spin_system,rho,npops);

% Report the state vector norm
report(spin_system,['state vector 2-norm: ' num2str(norm(rho,2))]);

% Locate npops most populated states and sort by amplitude
[~,sorting_index]=sort(abs(rho),1,'descend');
largest_elemts=rho(sorting_index(1:npops));

% Print the states and their populations
report(spin_system,[num2str(npops) ' most populated basis states (state, coeff, number)']);
for n=1:npops

    % Locate the hosting substance and its local descriptor row
    subst=find(sorting_index(n)<=spin_system.bas.offsets(2:end),1);
    row=sorting_index(n)-spin_system.bas.offsets(subst);
    spins=spin_system.chem.parts{subst};
    state_string=cell(1,numel(spins));
    for k=1:numel(spins)
        [L,M]=lin2lm(spin_system.bas.basis{subst}(row,k));
        if L==0
            state_string{k}='  ....  ';
        else
            state_string{k}=[' (' num2str(L,'%d') ',' num2str(M,'%+d') ') '];
        end
    end
    report(spin_system,['substance ' num2str(subst) ' spins [' num2str(spins) '] ' ...
                        cell2mat(state_string)              '    ' ...
                        num2str(largest_elemts(n),'%+5.3e') '    ' ...
                        num2str(sorting_index(n))]);
end

end

% Consistency enforcement
function grumble(spin_system,rho,npops)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('state analysis is only available for sphten-liouv formalism.');
end
if (~isnumeric(rho))||(~iscolumn(rho))
    error('rho parameter must be a column vector.');
end
if (~isnumeric(npops))||(~isscalar(npops))||(~isreal(npops))||...
   (npops<1)||(mod(npops,1)~=0)
    error('pops parameter must be a positive real integer.');
end
if npops>size(rho,1)
    error('number of populations exceeds state space dimension.');
end
end

% "One man with a dream will beat a hundred men with a grudge."
%
% A Russian saying

