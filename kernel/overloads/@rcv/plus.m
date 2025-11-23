% Adds RCV objects, sparse matrices or numeric scalars. Syntax:
%
%                         result=plus(a,b)
%
% Parameters:
%
%    a     - left operand (RCV, sparse or numeric scalar)
%    b     - right operand (RCV, sparse or numeric scalar)
%
% Outputs:
%
%    result- sum a+b
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/plus.m>

function result=plus(a,b)

% Check consistency
grumble(a,b);

% Add a scalar to an RCV object
if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)
    a.val=a.val+b;
    result=a;

elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)
    b.val=b.val+a;
    result=b;

% Add two RCV objects
elseif isa(a,'rcv')&&isa(b,'rcv')
    if a.numRows~=b.numRows||a.numCols~=b.numCols
        error('matrix sizes must match for addition.');
    end
    if a.isGPU&&~b.isGPU
        b.col=gpuArray(b.col);
        b.row=gpuArray(b.row);
        b.val=gpuArray(b.val);
    elseif ~a.isGPU&&b.isGPU
        a.col=gpuArray(a.col);
        a.row=gpuArray(a.row);
        a.val=gpuArray(a.val);
        a.isGPU=true;
    end
    a.row=[a.row;b.row];
    a.col=[a.col;b.col];
    a.val=[a.val;b.val];
    result=a;

% Add an RCV object and a sparse matrix
elseif isa(a,'rcv')&&issparse(b)
    if a.numRows~=size(b,1)||a.numCols~=size(b,2)
        error('matrix sizes must match for addition.');
    end
    b=rcv(b);
    result=plus(a,b);

elseif isa(b,'rcv')&&issparse(a)
    if b.numRows~=size(a,1)||b.numCols~=size(a,2)
        error('matrix sizes must match for addition.');
    end
    a=rcv(a);
    result=plus(a,b);
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

% Adam Smith was a genius, yes, but let's not
% forget he lived with his mum.
%
% Rory Sutherland
