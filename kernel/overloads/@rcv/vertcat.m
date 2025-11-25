% Vertical concatenation for RCV sparse matrices. Syntax:
%
%                      A=vertcat(A,B)
%
% Parameters:
%
%    A  - top RCV matrix
%    B  - bottom RCV matrix
%
% Outputs:
%
%    A  - concatenated RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/vertcat.m>

function A=vertcat(A,B)

% Check consistency
grumble(A,B);

% Align locations
if A.isGPU||B.isGPU
    A=gpuArray(A);
    B=gpuArray(B);
end

% Shift row indices
B.row=B.row+A.numRows;

% Concatenate indices
A.row=[A.row; B.row];
A.col=[A.col; B.col];
A.val=[A.val; B.val];

% Update row count in the result
A.numRows=A.numRows+B.numRows;

end

% Consistency enforcement
function grumble(A,B)
if ~isa(A,'rcv')||~isa(B,'rcv')
    error('both inputs must be RCV sparse matrices.');
end
if A.numCols~=B.numCols
    error('column counts must match for vertical concatenation.');
end
end

% Frankly speaking, my dear Karl, I do not like this modern word, which all
% weaklings use to cloak their feelings when they quarrel with the world
% because they do not possess, without labour or trouble, well-furnished
% palaces with vast sums of money and elegant carriages. This embitterment
% disgusts me and you are the last person from whom I would expect it. What
% grounds can you have for it? Has not everything smiled on you ever since
% your cradle? Has not nature endowed you with magnificent talents? Have
% not your parents lavished affection on you? Have you ever up to now been
% unable to satisfy your reasonable wishes? And have you not carried away
% in the most incomprehensible fashion the heart of a girl whom thousands
% envy you? Yet the first untoward event, the first disappointed wish,
% evokes embitterment! Is that strength? Is that a manly character?
%
% A letter to Karl Marx by his father in Nov 1837,
% Marx Engels Collected Works Vol 1, pp. 683-685.

