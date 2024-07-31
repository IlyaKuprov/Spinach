% Computes matrix exponential integrals of the following general type:
%
%             Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]
%
% Matrix A must be Hermitian. For further info see the paper by Char-
% les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
%
%                     R=expmint(spin_system,A,B,C,T)
%
% Parameters:
%
%    A,B,C   - the three matrices involved in the integral
%
%        T   - integration time
%
% Output:
%
%        R   - the resulting integral
%
% Note: the auxiliary matrix method is massively faster than either
%       commutator series or diagonalisation.
%
% Note: this is the most memory-intensive stage in a lot of calcula-
%       tions; memory recycling is aggressive.
% 
% ledwards@cbs.mpg.de
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=expmint.m>

function R=expmint(spin_system,A,B,C,T)

% Check consistency
grumble(A,B,C,T);

% Zero argument shortcut
if (T==0)||(nnz(B)==0), R=spalloc(size(B,1),size(B,2),0); return; end

% Build auxiliary matrix
auxmat=[-A, 1i*B; 0*A, -C];

% Get block extraction multipliers
BE1=spdiags(ones(size(A,1),1), 0        ,2*size(A,1),size(A,2));
BE2=spdiags(ones(size(A,1),1),-size(A,1),2*size(A,1),size(A,2));

% Reclaim memory
clear('A','B','C');

% Exponentiate the auxiliary matrix
auxmat=propagator(spin_system,auxmat,T);

% Inform the user
report(spin_system,'processing auxiliary matrix blocks...');

% Cut the blocks out
P=(BE1'*auxmat*BE1)'; Q=(BE1'*auxmat*BE2);

% Reclaim memory
clear('auxmat','BE1','BE2');

% Multiply up the blocks
R=clean_up(spin_system,P*Q,spin_system.tols.prop_chop);

% Reclaim memory
clear('P','Q');

end

% Consistency enforcement
function grumble(A,B,C,T)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||...
   (~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))
    error('A, B and C arguments must be matrices.');
end
if (~all(size(A)==size(B)))||(~all(size(B)==size(C)))
    error('A, B and C matrices must have the same dimension.');
end
if ~ishermitian(A)
    error('A matrix must be Hermitian.');
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% LADY NANCY ASTOR: "If you were my husband, Winston, I'd 
%                    put poison in your tea."
% WINSTON CHURCHILL: "If I were your husband, Nancy, I'd 
%                     drink it."

