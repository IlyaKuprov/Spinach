% Vertical concatenation for RCV sparse matrices. Syntax:
%
%                      A=vertcat(A,B,...)
%
% Parameters:
%
%    A,B,...   - RCV sparse matrices, top to bottom
%
% Outputs:
%
%    A         - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/vertcat.m>

function A=vertcat(varargin)

% Check consistency
grumble(varargin{:});

% Move all operands to the GPU if any of them is there
if any(cellfun(@(x)x.isGPU,varargin))
    varargin=cellfun(@gpuArray,varargin,'UniformOutput',false);
end

% Shift row indices by the running row count
rows=cell(nargin,1); cols=cell(nargin,1); vals=cell(nargin,1); nrows=int64(0);
for n=1:nargin
    rows{n}=varargin{n}.row+nrows;
    cols{n}=varargin{n}.col;
    vals{n}=varargin{n}.val;
    nrows=nrows+varargin{n}.numRows;
end

% Concatenate RCV arrays once
A=varargin{1};
A.row=vertcat(rows{:});
A.col=vertcat(cols{:});
A.val=vertcat(vals{:});
A.numRows=nrows;

end

% Consistency enforcement
function grumble(varargin)
if ~all(cellfun(@(x)isa(x,'rcv'),varargin))
    error('all inputs must be RCV sparse matrices.');
end
if numel(unique(cellfun(@(x)x.numCols,varargin)))>1
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

