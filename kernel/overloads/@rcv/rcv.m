% Creates an object of an RCV class. Syntax:
%
%                       obj=rcv(Mat)
%                       obj=rcv(dim1,dim2)
%                       obj=rcv(R,C,V)
%
% rcv (row, column, value) is a special format for storing and adding large
% sparse matrices.
%
% Parameters:
%
%    Mat    - sparse or numeric matrix
%    dim1   - number of rows for an empty object
%    dim2   - number of columns for an empty object
%    R      - row indices of non-zero entries
%    C      - column indices of non-zero entries
%    V      - values corresponding to entries in R and C
%
% Outputs:
%
%    obj    - an object that behaves like a sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/rcv.m>

classdef (InferiorClasses={?gpuArray}) rcv

    properties
        row=int64([]);     % Row indices (int64)
        col=int64([]);     % Column indices (int64)
        val=double([]);    % Values (double)
        numRows=int64(0);  % Number of rows in the original matrix
        numCols=int64(0);  % Number of columns in the original matrix
        isGPU=false;       % Flag indicating if stored on GPU
    end

    methods

        function obj=rcv(varargin)

            % Check consistency
            grumble_constructor(varargin{:});

            if nargin==1

                % Return the input if it is already RCV
                input=varargin{1};
                if isa(input,'rcv')
                    obj=input;

                % Convert numeric or sparse input to RCV
                elseif isnumeric(input)||issparse(input)||isa(input,'gpuArray')
                    obj.isGPU=isa(input,'gpuArray');
                    obj.numRows=int64(size(input,1));
                    obj.numCols=int64(size(input,2));
                    [row,col,val]=find(input);
                    if obj.isGPU
                        obj.row=gpuArray(int64(row(:)));
                        obj.col=gpuArray(int64(col(:)));
                        obj.val=gpuArray(double(val(:)));
                    else
                        obj.row=int64(row(:));
                        obj.col=int64(col(:));
                        obj.val=double(val(:));
                    end
                end

            elseif nargin==2

                % Create an empty object with specified dimensions
                obj.numRows=int64(varargin{1});
                obj.numCols=int64(varargin{2});
                obj.row=int64([]);
                obj.col=int64([]);
                obj.val=double([]);
                obj.isGPU=false;

            elseif nargin==3

                % Build an object from row, column and value vectors
                obj.row=int64(varargin{1}(:));
                obj.col=int64(varargin{2}(:));
                obj.val=double(varargin{3}(:));
                obj.numRows=max(obj.row);
                obj.numCols=max(obj.col);
                obj.isGPU=isa(varargin{1},'gpuArray')||isa(varargin{2},'gpuArray')||isa(varargin{3},'gpuArray');
                if obj.isGPU
                    obj.row=gpuArray(obj.row);
                    obj.col=gpuArray(obj.col);
                    obj.val=gpuArray(obj.val);
                end
            end
        end

        function answer=isnumeric(obj) %#ok<MANU>

            % Check consistency
            grumble_unary(obj);

            % RCV matrices are numeric
            answer=true();

        end

        function answer=ismatrix(obj) %#ok<MANU>

            % Check consistency
            grumble_unary(obj);

            % RCV matrices are matrices
            answer=true();

        end

        function answer=isfloat(obj) %#ok<MANU>

            % Check consistency
            grumble_unary(obj);

            % RCV matrices are floats
            answer=true();

        end

    end

end

% Consistency enforcement for constructor
function grumble_constructor(varargin)
if (nargin<1)||(nargin>3)
    error('use one, two or three arguments to create an rcv object.');
end
if nargin==1
    input=varargin{1};
    if ~(isa(input,'rcv')||isnumeric(input)||issparse(input)||isa(input,'gpuArray'))
        error('single argument must be an rcv, numeric or sparse matrix.');
    end
elseif nargin==2
    dim1=varargin{1}; dim2=varargin{2};
    if ~(isnumeric(dim1)&&isnumeric(dim2)&&isscalar(dim1)&&isscalar(dim2))
        error('two arguments must be numeric scalars specifying dimensions.');
    end
elseif nargin==3
    row=varargin{1}; col=varargin{2}; val=varargin{3};
    if ~(isnumeric(row)&&isnumeric(col)&&isnumeric(val))
        error('row, column and value arrays must be numeric.');
    end
    if ~(isvector(row)&&isvector(col)&&isvector(val))
        error('row, column and value inputs must be vectors.');
    end
    if ~(numel(row)==numel(col)&&numel(col)==numel(val))
        error('row, column and value vectors must be the same length.');
    end
end
end

% Consistency enforcement for unary checks
function grumble_unary(obj)
if ~isa(obj,'rcv')
    error('the input must be an rcv object.');
end
if ~isscalar(obj)
    error('the input must be a scalar rcv object.');
end
end

% There are weeks when decades happen.
%
% Vladimir Lenin
