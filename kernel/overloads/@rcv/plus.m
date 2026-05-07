% Adds things to RCV sparse matrices. Syntax:
%
%                  C=plus(A,B)
%
% Parameters:
%
%    A  - left operand
%
%    B  - right operand
%
% Outputs:
%
%    C  - sum A+B, RCV sparse matrix
%
% m.keitel@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rcv/plus.m>

function C=plus(A,B)

% Check consistency
grumble(A,B);

% Process the special cases
if isa(A,'rcv')&&isscalar(B)&&isnumeric(B)

    % Explain the refusal to add a scalar to an RCV object
    error('adding a scalar makes a matrix non-sparse; inflate manually first.');

% Add a scalar to an RCV sparse matrix
elseif isa(B,'rcv')&&isscalar(A)&&isnumeric(A)
    
    % Explain the refusal to add a scalar to an RCV object
    error('adding a scalar makes a matrix non-sparse; inflate manually first.');

% Add two RCV sparse matrices
elseif isa(A,'rcv')&&isa(B,'rcv')

    % Check for dimension match
    if A.numRows~=B.numRows||A.numCols~=B.numCols
        error('matrix sizes must match.');
    end

    % Align locations
    if A.isGPU||B.isGPU
        A=gpuArray(A);
        B=gpuArray(B);
    end
    
    % RCV addition
    A.row=[A.row; B.row];
    A.col=[A.col; B.col];
    A.val=[A.val; B.val]; C=A;

% RCV sparse + Matlab sparse
elseif isa(A,'rcv')&&issparse(B)

    % Check for dimension match
    if A.numRows~=size(B,1)||A.numCols~=size(B,2)
        error('matrix sizes must match.');
    end

    % Recursive call
    C=plus(A,rcv(B));

elseif isa(B,'rcv')&&issparse(A)

    % Check for dimension match
    if B.numRows~=size(A,1)||B.numCols~=size(A,2)
        error('matrix sizes must match.');
    end

    % Recursive call
    C=plus(rcv(A),B);

end

end

% Consistency enforcement
function grumble(A,B)
if (~isa(A,'rcv'))&&(~isa(B,'rcv'))
    error('at least one input must be an RCV sparse matrix.');
end
if isa(A,'rcv')&&(~isa(B,'rcv'))&&(~issparse(B))&&((~isnumeric(B))||(~isscalar(B)))
    error('the second input must be an RCV sparse matrix, a Matlab sparse matrix, or a numeric scalar.');
end
if isa(B,'rcv')&&(~isa(A,'rcv'))&&(~issparse(A))&&((~isnumeric(A))||(~isscalar(A)))
    error('the first input must be an RCV sparse matrix, a Matlab sparse matrix, or a numeric scalar.');
end
end

% Adam Smith was a genius, yes, but let's not
% forget he lived with his mum.
%
% Rory Sutherland

