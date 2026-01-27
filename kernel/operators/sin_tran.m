% Single transition operators, spanning the space of matri-
% ces of the specified dimension. The set is returned as a
% cell array of sparse matrices. The first element is the
% unit matrix, the rest use serpentine indexing where the
% position in the cell array maps in the following way to
% the location of a single non-zero in the matrix:
%
%                  (2)  (4)   (7)   (11) 
%                  (3)  (6)   (10)  (14) 
%                  (5)  (9)   (13)  (16)
%                  (8)  (12)  (15)   xx
%
% and likewise for larger matrices. The last element is ab-
% sent because it makes the set linearly dependent. Syntax:
%
%                     A=sin_tran(dim)
%
% Parameters:
%
%   dim - dimension of the matrices
%   
% Outputs:
%
%     A - a cell array of matrices, structured
%         as described above; matrices are re-
%         turned as complex to avoid expensive
%         reallocations later
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sin_tran.m>

function A=sin_tran(dim)

% Check consistency
grumble(dim);

% Empty array
A=cell(dim^2,1);

% First element is unit
A{1}=complex(speye(dim));

% Fill the array
parfor n=1:(dim^2-1)
    [k,q]=lin2kq(dim,n-1);
    A{n+1}=spalloc(dim,dim,1);
    A{n+1}(k+1,q+1)=1;
    A{n+1}=complex(A{n+1});
end

end

% Consistency enforcement
function grumble(dim)
if (~isnumeric(dim))||(~isscalar(dim))||...
   (~isreal(dim))||(dim<2)||(mod(dim,1)~=0)
    error('dim must be a real integer greater than 1.');
end
end

% Have nothing in your house that you do not know
% to be useful or believe to be beautiful.
%
% William Morris

