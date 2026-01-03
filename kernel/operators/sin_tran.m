% Single transition operators, numbered in the following
% order by where a single unit element appears:
%
%                       |1 2 5 ...|
%                       |4 3 6 ...|
%                       |9 8 7 ...|
%                       |.........|
% Syntax:
%
%                     A=sin_tran(dim,n)
%
% Parameters:
%
%   dim -  the dimension of the resulting
%          sparse square matrix
%   
%   n   -  the position of the unit element
%          by the above numbering scheme, 0
%          produces a unit matrix
%
% Outputs:
%
%   A   -  a sparse matrix, complex double
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sin_tran.m>

function A=sin_tran(dim,n)

% Check consistency
grumble(dim,n);

% Unit matrix shortcut
if n==0

    % Complex unit matrix
    A=complex(speye(dim));

else

    % Do the clockwise spiral indexing
    k=ceil(sqrt(double(n))); t=n-(k-1)^2;
    row=min(t,k); col=k-max(t-k,0);

    % Build the complex matrix
    A=sparse(row,col,1,dim,dim);
    A=complex(A);

end

end

% Consistency enforcement
function grumble(dim,n)
if (~isnumeric(dim))||(~isscalar(dim))||...
   (~isreal(dim))||(dim<1)||(mod(dim,1)~=0)
    error('dim must be a positive real integer.');
end
if (~isnumeric(n))||(~isscalar(n))||...
   (~isreal(n))||(n<0)||(mod(n,1)~=0)
    error('dim must be a non-negative real integer.');
end
if n>dim^2
    error('n cannot exceed dim^2');
end
end

% Have nothing in your house that you do not know
% to be useful or believe to be beautiful.
%
% William Morris

