% Adds things to RCV sparse matrices. Syntax:
%
%                result=plus(a,b)
%
% Parameters:
%
%    a  - left operand
%
%    b  - right operand
%
% Outputs:
%
%    result - sum a+b
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/plus.m>

function result=plus(a,b)

% Check consistency
grumble(a,b);

% Add a scalar to an RCV sparse matrix
if isa(a,'rcv')&&isscalar(b)&&isnumeric(b)

    a.val=a.val+b; result=a;

% Add a scalar to an RCV sparse matrix
elseif isa(b,'rcv')&&isscalar(a)&&isnumeric(a)
    
    b.val=b.val+a; result=b;

% Add two RCV sparse matrices
elseif isa(a,'rcv')&&isa(b,'rcv')

    % Check for dimension match
    if a.numRows~=b.numRows||a.numCols~=b.numCols
        error('matrix sizes must match for addition.');
    end

    % Location matching
    if a.isGPU&&(~b.isGPU)
        
        % GPU upload
        b.row=gpuArray(b.row);
        b.col=gpuArray(b.col);
        b.val=gpuArray(b.val);
        b.isGPU=true;

    elseif (~a.isGPU)&&b.isGPU

        % GPU upload
        a.col=gpuArray(a.col);
        a.row=gpuArray(a.row);
        a.val=gpuArray(a.val);
        a.isGPU=true;

    end
    
    % RCV addition
    a.row=[a.row; b.row];
    a.col=[a.col; b.col];
    a.val=[a.val; b.val];
    result=a;

% RCV sparse + Matlab sparse
elseif isa(a,'rcv')&&issparse(b)

    % Check for dimension match
    if a.numRows~=size(b,1)||a.numCols~=size(b,2)
        error('matrix sizes must match for addition.');
    end

    % Recursive call
    result=plus(a,rcv(b));

elseif isa(b,'rcv')&&issparse(a)

    % Check for dimension match
    if b.numRows~=size(a,1)||b.numCols~=size(a,2)
        error('matrix sizes must match for addition.');
    end

    % Recursive call
    result=plus(rcv(a),b);

end

end

% Consistency enforcement
function grumble(a,b)

% Needs some thinking

end

% Adam Smith was a genius, yes, but let's not
% forget he lived with his mum.
%
% Rory Sutherland

