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
%   in the basis set.
%
% Note: this function requires a spherical tensor basis set.
%
% i.kuprov@soton.ac.uk
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
largest_states=spin_system.bas.basis(sorting_index(1:npops),:);

% Print the states and their populations
report(spin_system,[num2str(npops) ' most populated basis states (state, coeff, number)']);
for n=1:npops
    state_string=cell(1,spin_system.comp.nspins);
    for k=1:spin_system.comp.nspins
        [l,m]=lin2lm(largest_states(n,k));
        if l==0
            state_string{k}='  ....  ';
        else
            state_string{k}=[' (' num2str(l,'%d') ',' num2str(m,'%+d') ') '];
        end
    end
    disp([cell2mat(state_string) '    ' num2str(largest_elemts(n),'%+5.3e') '    ' num2str(sorting_index(n))]);
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

