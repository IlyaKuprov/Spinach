% Multiplication for RCV sparse matrices. Syntax:
%
%                     C=mtimes(A,B)
%
% Parameters:
%
%    A    - left operand 
%
%    B    - right operand 
%
% Outputs:
%
%    C    - product A*B as a Matlab sparse matrix if
%           both operands are RCV or Matlab matrices
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/mtimes.m>

function C=mtimes(A,B)

% Check consistency
grumble(A,B);

% RCV sparse by a scalar
if isa(A,'rcv')&&isscalar(B)&&isnumeric(B)

    A.val=A.val*B; C=A;

% Scalar by RCV sparse
elseif isa(B,'rcv')&&isscalar(A)&&isnumeric(A)

    B.val=B.val*A; C=B;

% RCV sparse by RCV sparse
elseif isa(A,'rcv')&&isa(B,'rcv')

    % Check dimension consistency
    if A.numCols~=B.numRows
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    C=sparse(A)*sparse(B);

elseif isa(A,'rcv')&&issparse(B)

    % Check dimension consistency
    if A.numCols~=size(B,1)
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    C=sparse(A)*B;

% Matlab sparse by RCV sparse
elseif issparse(A)&&isa(B,'rcv')

    % Check dimension consistency
    if size(A,2)~=B.numRows
        error('inner matrix dimensions must agree.');
    end

    % Result is Matlab sparse
    C=A*sparse(B);

end

end

% Consistency enforcement
function grumble(A,B)

% Needs some thought

end

% When the Earl of Sandwich, speaking in parliament, told John Wilkes that
% the latter would either die on the gallows or of the pox, Wilkes politely
% responded that it would depend on whether he embraced the earl's princi-
% ples or his mistress.
%
% Taki Theodoracopoulos

