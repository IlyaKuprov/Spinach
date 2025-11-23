% Creates an object of a RCV class. Syntax:
%
%                          p=rcv(Mat) or
%                          p=rcv(Dim1,Dim2) or
%                          p=rcv(R,C,V)
%
% rcv (row, column, value) is a special format for storing and adding large
% sparse matrices.
%
% Parameters:
% 
%        Mat - Sparse Matrix
%        Dim1,Dim2 - Empty Matrix Object with Size Dim1xDim2
%        R,C,V - Matrix with entries in rows R, columns C, with values V
%
% Outputs:
%
%        rcv - An object that behaves like a sparse matrix
%
% m.keitel@soton.ac.uk

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

            if nargin==1

                % Case: Single argument
                input=varargin{1};
                if isa(input,'rcv')
                    obj=input;
                elseif isnumeric(input)||issparse(input)||isa(input,'gpuArray')
                    % Convert matrix to RCV format
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
                else
                    error('RCV:InvalidInput','Single argument must be a numeric matrix.');
                end

            elseif nargin==2

                % Case: Two arguments (numRows, numCols)
                numRows=varargin{1};
                numCols=varargin{2};
                if isnumeric(numRows)&&isnumeric(numCols)&&...
                   isscalar(numRows)&&isscalar(numCols)
                    obj.numRows=int64(numRows);
                    obj.numCols=int64(numCols);
                    obj.row=int64([]);
                    obj.col=int64([]);
                    obj.val=double([]);
                    obj.isGPU=false;
                else
                    error('RCV:InvalidInput','Two arguments must be numeric scalars specifying matrix dimensions.');
                end

            elseif nargin==3
                % Case: Three arguments (col, row, val vectors)
                row=varargin{1};
                col=varargin{2};
                val=varargin{3};

                if isnumeric(row)&&isnumeric(col)&&isnumeric(val)&&...
                   isvector(row)&&isvector(col)&&isvector(val)&&...
                   numel(row)==numel(col)&&numel(col)==numel(val)

                    obj.row=int64(row(:));
                    obj.col=int64(col(:));
                    obj.val=double(val(:));
                    obj.numRows=max(int64(row)); % Infer dimensions from max indices
                    obj.numCols=max(int64(col));
                    obj.isGPU=false; % Assume CPU unless vectors are GPU arrays

                    if isa(col,'gpuArray')||isa(row,'gpuArray')||isa(val,'gpuArray')
                        obj.row=gpuArray(obj.row);
                        obj.col=gpuArray(obj.col);
                        obj.val=gpuArray(obj.val);
                        obj.isGPU=true;
                    end
                else
                    error('RCV:InvalidInput','Three arguments must be numeric vectors of the same length.');
                end

            else
                error('RCV:InvalidInput','Invalid number of arguments. Use 1 (matrix), 2 (dimensions), or 3 (CRV vectors).');
            end
        end

        % RCV matrices are numeric
        function answer=isnumeric(obj) %#ok<MANU>
    
            answer=true(); % Always
    
        end
    
        % RCV matrices are matrices
        function answer=ismatrix(obj) %#ok<MANU>
    
            answer=true(); % Always
    
        end
    
        % RCV matrices are floats
        function answer=isfloat(obj) %#ok<MANU>
    
            answer=true(); % Always
    
        end

    end

end

% There are weeks when decades happen.
%
% Vladimir Lenin

