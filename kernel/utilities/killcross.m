% Zeroes the specified rows and columns of a matrix. Syntax:
%
%                  M=killcross(M,f1idx,f2idx)
%
% Parameters:
%
%     M      - a matrix
%
%     f1idx  - numbers of the columns that 
%              should be zeroed
%
%     f2idx  - numbers of the rows that 
%              should be zeroed
%
% Outputs:
%
%     M      - a matrix
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=killcross.m>

function M=killcross(M,f1idx,f2idx)

% Check consistency
grumble(M,f1idx,f2idx);

% Wipe the indices
M(f2idx,:)=0; M(:,f1idx)=0;

end

% Consistency enforcement
function grumble(M,f1idx,f2idx)
if (~isnumeric(M))||(~ismatrix(M))
    error('M must be a matrix.');
end
if (~isnumeric(f1idx))||(~isnumeric(f2idx))||...
   (~isreal(f1idx))||(~isreal(f2idx))||...
   any(f1idx<1)||any(f2idx<1)||any(mod(f1idx,1)~=0)||...
   any(mod(f2idx,1)~=0)
    error('index arrays must contain positive integers.');
end
if (~isequal(f1idx,unique(f1idx)))||(~isequal(f2idx,unique(f2idx)))
    error('repeated elements not allowed in the index arrays.');
end
if any(f1idx>size(M,2))||any(f2idx>size(M,1))
    error('index array element exceeds spectrum matrix dimension.');
end
end

% A narcissist is someone better-looking than you are.
%
% Gore Vidal

