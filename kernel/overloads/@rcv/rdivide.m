% Divides an RCV sparse matrix by a numeric scalar. Syntax:
%
%                       A=rdivide(A,k)
%
% Parameters:
%
%    A    - RCV sparse matrix
%
%    k - numeric scalar
%
% Outputs:
%
%    A    - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/rdivide.m>

function A=rdivide(A,k)

% Check consistency
grumble(A,k);

% Divide stored values by the scalar
A.val=A.val/k;

end

% Consistency enforcement
function grumble(A,k)
if ~isa(A,'rcv')
    error('the first argument must be an RCV sparse matrix.');
end
if (~isnumeric(k))||(~isscalar(k))
    error('division is only defined for numeric scalars.');
end
end

% Главное памятник поставить, а голуби сами прилетят.
%
% Святослав Вернидубович Кривич

