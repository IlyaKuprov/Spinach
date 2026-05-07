% Sparse block diagonal matrix from a stack of matrix blocks.
% Syntax:
%
%                 S=sp_block_diag(A)
%                 S=sp_block_diag(A,B,C,...)
%
% Parameters:
%
%    A       - a floating point array; in the first syntax, the
%              first two dimensions are matrix dimensions and the
%              remaining dimensions enumerate the blocks
%
%    A,B,C   - floating point matrices to be placed on the block
%              diagonal in the second syntax
%
% Outputs:
%
%    S       - sparse block diagonal matrix
%
% Notes: this function is a Spinach-local replacement for Matlab's
%        spblkdiag function from the Model-Based Calibration toolbox.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sp_block_diag.m>

function S=sp_block_diag(varargin)

% Check consistency
if nargin==0
    error('at least one input array must be supplied.');
end
grumble(varargin{:});

% Single input: block stack
if nargin==1

    % Get the array
    A=varargin{1};

    % Get block dimensions
    nrows=size(A,1); ncols=size(A,2);
    block_dims=size(A);
    if numel(block_dims)<3
        nblocks=1;
    else
        nblocks=prod(block_dims(3:end));
    end

    % Reshape higher dimensions into a block index
    A=reshape(A,[nrows ncols nblocks]);

    % Get non-zero elements
    [lin_idx,~,vals]=find(A(:));
    [rows,cols,blocks]=ind2sub([nrows ncols nblocks],lin_idx);

    % Assemble the matrix
    S=sparse(rows+(blocks-1)*nrows,cols+(blocks-1)*ncols,vals,...
             nrows*nblocks,ncols*nblocks);

else

    % Count block dimensions
    nrows=zeros(nargin,1); ncols=zeros(nargin,1);
    for n=1:nargin
        nrows(n)=size(varargin{n},1);
        ncols(n)=size(varargin{n},2);
    end

    % Compute offsets
    row_offsets=[0; cumsum(nrows(1:(end-1)))];
    col_offsets=[0; cumsum(ncols(1:(end-1)))];

    % Preallocate index cells
    rows=cell(nargin,1); cols=cell(nargin,1); vals=cell(nargin,1);

    % Loop over blocks
    for n=1:nargin
        [rows{n},cols{n},vals{n}]=find(varargin{n});
        rows{n}=rows{n}(:)+row_offsets(n);
        cols{n}=cols{n}(:)+col_offsets(n);
        vals{n}=vals{n}(:);
    end

    % Assemble the matrix
    S=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),...
             sum(nrows),sum(ncols));

end

end

% Consistency enforcement
function grumble(varargin)
for n=1:nargin
    if (~isfloat(varargin{n}))||(~isnumeric(varargin{n}))
        error('all inputs must be double or single arrays.');
    end
    if (nargin>1)&&(~ismatrix(varargin{n}))
        error('with multiple inputs, all inputs must be matrices.');
    end
end
end

% Hopeless? It's not hopeless.
% Doubtful, but not hopeless.
%
% A-ha

