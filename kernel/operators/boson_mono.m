% Bosonic monomial operators of the following structure:
%
%                  B(k,q)=(Cr^k)*(An^q)
%
% obeying the following commutation relations with the po-
% pulation number operator N:
%
%                [N,B(k,q)]=(k-q)*B(k,q)
% 
% Syntax:
%
%                 B=boson_mono(nlevels)
%
% Parameters:
%
%     nlevels - number of bosonic ladder population 
%               levels, k and q go from 0 to nlevels-1
%
% Outputs:
%
%       B - a cell array with the following numbering
%           map between (k,q) and a single index:
%
%                (0,0)(0,1)(0,2)     (1)(3)(6)
%                (1,0)(1,1)(1,2) <=> (2)(5)(8)
%                (2,0)(2,1)(2,2)     (4)(7)(9)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=boson_mono.m>

function B=boson_mono(nlevels)

% Check consistency
grumble(nlevels);

% Get generators
A=weyl(nlevels);

% Generate monomials 
B=cell(nlevels,nlevels);
for k=0:(nlevels-1)
    for q=0:(nlevels-1)
        B{k+1,q+1}=(A.c^k)*(A.a^q);
    end
end

% Build the serpentine index
[rows,cols]=ndgrid(1:nlevels);
[~,idx]=sortrows([rows(:)+cols(:), -rows(:)]);

% Arrange
B=B(idx);

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('nlevels must be a positive integer.');
end
end


% "Oh, so you really exist. I thought Littlewood 
%  was a pseudonym Hardy used for his less impor-
%  tant articles."
%
% Norbert Wiener, to John Littlewood

