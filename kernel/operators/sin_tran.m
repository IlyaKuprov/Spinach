% Single transition operators, spanning the space of matri-
% ces of the specified dimension. The set is returned as a
% cell array of sparse matrices using serpentine indexing
% where the position in the cell array maps in the follow-
% ing way to the location of the single non-zero:
%
%                  (1)  (3)   (6)   (10) 
%                  (2)  (5)   (9)   (13) 
%                  (4)  (8)   (12)  (15)
%                  (7)  (11)  (14)  (16)
%
% and likewise for larger matrices. Syntax:
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

% Fill the array
parfor n=1:dim^2

    % Unit element indices
    [k,q]=lin2kq(dim,n);

    % Matrix construction
    A{n}=sparse(k,q,1,dim,dim);

    % Complex type
    A{n}=complex(A{n});

end

end

% Consistency enforcement
function grumble(dim)
if (~isnumeric(dim))||(~isscalar(dim))||...
   (~isreal(dim))||(dim<1)||(mod(dim,1)~=0)
    error('dim must be a positive real integer.');
end
end

% Have nothing in your house that you do not know
% to be useful or believe to be beautiful.
%
% William Morris

