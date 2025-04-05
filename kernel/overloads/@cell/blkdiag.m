% Block-diagonal cell array from two cell arrays, all
% other elements are set to empty cells. Syntax:
%
%                   C=blkdiag(A,B)
%
% Parameters:
%
%    A,B - cell arrays
%
% Outputs:
%
%    C   - cell array
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=blkdiag.m>

function C=blkdiag(A,B)

% Decide the dimensions
dim_a=size(A); dim_b=size(B);

% Make an empty array
C=cell(dim_a+dim_b);

% Fill in the blocks
C(1:dim_a(1),1:dim_a(2))=A;
C((dim_a(1)+1):end,(dim_a(2)+1):end)=B;

end

% Мы за мир, но есть нюансы.
%
% Владимир Путин

