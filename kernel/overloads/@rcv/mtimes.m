% Multiplication for RCV sparse matrices. Syntax:
%
%                       result=mtimes(a,b)
%
% Parameters:
%
%    a     - left operand 
%
%    b     - right operand 
%
% Outputs:
%
%    result - product a*b
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/mtimes.m>

function result=mtimes(a,b)

% Check consistency
grumble(a,b);

% RCV sparse by a scalar
if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)

    a.val=a.val*b; result=a;

% Scalar by RCV sparse
elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)

    b.val=b.val*a; result=b;

% RCV sparse by RCV sparse
elseif isa(a,'rcv')&&isa(b,'rcv')

    % Check dimension consistency
    if a.numCols~=b.numRows
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    result=sparse(a)*sparse(b);

elseif isa(a,'rcv')&&issparse(b)

    % Check dimension consistency
    if a.numCols~=size(b,1)
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    result=sparse(a)*b;

% Matlab sparse by RCV sparse
elseif issparse(a)&&isa(b,'rcv')

    % Check dimension consistency
    if size(a,2)~=b.numRows
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    result=a*sparse(b);

end

end

% Consistency enforcement
function grumble(a,b)

% Needs some thought

end

% When the Earl of Sandwich, speaking in parliament, told John Wilkes that
% the latter would either die on the gallows or of the pox, Wilkes politely
% responded that it would depend on whether he embraced the earl's princi-
% ples or his mistress.
%
% Taki Theodoracopoulos

