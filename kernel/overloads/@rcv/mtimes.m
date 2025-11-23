% Matrix multiplication for RCV objects. Syntax:
%
%                       result=mtimes(a,b)
%
% Parameters:
%
%    a     - left operand (RCV, sparse or scalar numeric)
%    b     - right operand (RCV, sparse or scalar numeric)
%
% Outputs:
%
%    result- product a*b
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/mtimes.m>

function result=mtimes(a,b)

% Check consistency
grumble(a,b);

% Scalar multiplication with an RCV object
if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)
    a.val=a.val*b;
    result=a;

elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)
    b.val=b.val*a;
    result=b;

% Multiplication between two RCV objects
elseif isa(a,'rcv')&&isa(b,'rcv')
    if a.numCols~=b.numRows
        error('inner matrix dimensions must agree.');
    end
    result=sparse(a)*sparse(b);

% Multiplication between RCV and sparse matrices
elseif isa(a,'rcv')&&issparse(b)
    if a.numCols~=size(b,1)
        error('inner matrix dimensions must agree.');
    end
    result=sparse(a)*b;

elseif issparse(a)&&isa(b,'rcv')
    if size(a,2)~=b.numRows
        error('inner matrix dimensions must agree.');
    end
    result=a*sparse(b);
end

end

% Consistency enforcement
function grumble(a,b)
if ~(isa(a,'rcv')||isa(b,'rcv')||issparse(a)||issparse(b)||isnumeric(a)||isnumeric(b))
    error('inputs must be rcv, sparse or numeric objects.');
end
if isa(a,'rcv')&&~isscalar(a)
    error('rcv objects must be scalar.');
end
if isa(b,'rcv')&&~isscalar(b)
    error('rcv objects must be scalar.');
end
if isnumeric(a)&&~isscalar(a)
    error('numeric inputs must be scalar.');
end
if isnumeric(b)&&~isscalar(b)
    error('numeric inputs must be scalar.');
end
end

% When the Earl of Sandwich, speaking in parliament, told John Wilkes that
% the latter would either die on the gallows or of the pox, Wilkes politely
% responded that it would depend on whether he embraced the earl's princi-
% ples or his mistress.
%
% Taki Theodoracopoulos
