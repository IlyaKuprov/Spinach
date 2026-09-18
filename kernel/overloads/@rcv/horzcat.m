% Horizontal concatenation for RCV sparse matrices. Syntax:
%
%                        A=horzcat(A,B,...)
%
% Parameters:
%
%    A,B,...   - RCV sparse matrices, left to right
%
% Outputs:
%
%    A         - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/horzcat.m>

function A=horzcat(varargin)

% Check consistency
grumble(varargin{:});

% Move all operands to the GPU if any of them is there
if any(cellfun(@(x)x.isGPU,varargin))
    varargin=cellfun(@gpuArray,varargin,'UniformOutput',false);
end

% Shift column indices by the running column count
rows=cell(nargin,1); cols=cell(nargin,1); vals=cell(nargin,1); ncols=int64(0);
for n=1:nargin
    rows{n}=varargin{n}.row;
    cols{n}=varargin{n}.col+ncols;
    vals{n}=varargin{n}.val;
    ncols=ncols+varargin{n}.numCols;
end

% Concatenate RCV arrays once
A=varargin{1};
A.row=vertcat(rows{:});
A.col=vertcat(cols{:});
A.val=vertcat(vals{:});
A.numCols=ncols;

end

% Consistency enforcement
function grumble(varargin)
if ~all(cellfun(@(x)isa(x,'rcv'),varargin))
    error('all inputs must be RCV sparse matrices.');
end
if numel(unique(cellfun(@(x)x.numRows,varargin)))>1
    error('row counts must match for horizontal concatenation.');
end
end

% The back half of your forties is a cursed age. It's not
% so much that (as the cliché goes) the policemen look
% younger - it's more that you're gloomily aware that
% you haven't done anything that they'd even consider
% arresting you for.
%
% Sam Leith

