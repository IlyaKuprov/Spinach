% Divides an RCV sparse matrix by a numeric scalar. Syntax:
%
%                   obj=rdivide(obj,scalar)
%
% Parameters:
%
%    obj    - RCV sparse matrix
%
%    scalar - numeric scalar
%
% Outputs:
%
%    obj    - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/rdivide.m>

function obj=rdivide(obj,scalar)

% Check consistency
grumble(obj,scalar);

% Divide stored values by the scalar
obj.val=obj.val/scalar;

end

% Consistency enforcement
function grumble(obj,scalar)
if ~isa(obj,'rcv')
    error('the first argument must be an RCV sparse matrix.');
end
if (~isnumeric(scalar))||(~isscalar(scalar))
    error('division is only defined for numeric scalars.');
end
end

% Главное памятник поставить, а голуби сами прилетят.
%
% Святослав Вернидубович Кривич

