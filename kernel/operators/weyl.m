% Weyl boson operators (sparse, see below for normalisa-
% tion convention) for a bosonic mode with a user-speci-
% fied population number truncation. Syntax:
%
%                      A=weyl(nlevels)
%
% Parameters:
%
%     nlevels - an integer specifying the 
%               number of population levels
%
% Outputs:
%
%     A.u - unit operator
%
%     A.c - creation operator
%
%     A.a - annihilation operator
%
%     A.n - population number operator
%
% Note: the matrices are normalised to obey the following
%       relations for all energy level counts
%
%                  A.c*A.a=A.n,  [A.n,A.c]=A.c
%                [A.n,A.a]=-A.a,  [A.a,A.c]=A.u
%
%       except for the edge state at which [A.a,A.c] ele-
%       ment is (1-nlevels), this is unavoidable.
%
% Note: arrays are declared complex at build time to avoid 
%       expensive reallocation operations later on.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=weyl.m>

function A=weyl(nlevels)

% Ensure internal consistency
grumble(nlevels); nlevels=double(nlevels);

% Creation operator
diags=sqrt(1:nlevels); 
A.c=spdiags(diags',-1,nlevels,nlevels);
A.c=complex(A.c);

% Number operator
diags=0:(nlevels-1);
A.n=spdiags(diags',0,nlevels,nlevels);
A.n=complex(A.n);

% Annihilation operator
diags=sqrt(0:(nlevels-1));
A.a=spdiags(diags',+1,nlevels,nlevels);
A.a=complex(A.a);

% Unit operator
A.u=complex(speye(nlevels));
    
end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(mod(nlevels,1)~=0)||(nlevels<1)
    error('nlevels must be a positive real integer.');
end
end

% I happen to be a physicist who started life as a 
% mathematician. As a working physicist, I am acutely
% aware of the fact that the marriage between mathema-
% tics and physics, which was so enormously fruitful
% in past centuries, has recently ended in divorce.
%
% Freeman Dyson

