% Creates an RCV (row-column-value storage) sparse matrix. Syntax:
%
%                       obj=rcv(M)
%                       obj=rcv(dim1,dim2)
%                       obj=rcv(R,C,V,dim1,dim2)
%
% Parameters:
%
%    M      - a Matlab matrix
%
%    dim1   - number of rows 
%
%    dim2   - number of columns 
%
%    R      - row indices of non-zero entries
%
%    C      - column indices of non-zero entries
%
%    V      - values corresponding to entries in R and C
%
% Outputs:
%
%    obj    - an RCV sparse matrix object
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
            grumble(varargin{:});

            if nargin==1

                % Do nothing if already RCV
                if isa(varargin{1},'rcv')
                    
                    % Return input
                    obj=varargin{1};

                % Convert Matlab matrices to RCV
                elseif ismatrix(input)
                    
                    % Preserve location
                    obj.isGPU=isa(input,'gpuArray');

                    % Get dimensions
                    obj.numRows=int64(size(input,1));
                    obj.numCols=int64(size(input,2));

                    % Get non-zeroes
                    [row,col,val]=find(input);

                    % Make the object
                    obj.row=int64(row(:));
                    obj.col=int64(col(:));
                    obj.val=double(val(:));
                   
                else

                    % Complain and bomb out
                    error('input cannot be converted into RCV sparse format.');
                
                end

            elseif nargin==2

                % Empty object with specified dimensions
                obj.numRows=int64(varargin{1});
                obj.numCols=int64(varargin{2});
                obj.row=int64([]);
                obj.col=int64([]);
                obj.val=double([]);
                obj.isGPU=false;

            elseif nargin==5

                % Build an object from scratch
                obj.row=int64(varargin{1}(:));
                obj.col=int64(varargin{2}(:));
                obj.val=double(varargin{3}(:));
                obj.numRows=int64(varargin{4});
                obj.numCols=int64(varargin{5});

                % Decide the location
                obj.isGPU=isa(varargin{1},'gpuArray')||...
                          isa(varargin{2},'gpuArray')||...
                          isa(varargin{3},'gpuArray');

                % Upload to GPU
                if obj.isGPU
                    obj.row=gpuArray(obj.row);
                    obj.col=gpuArray(obj.col);
                    obj.val=gpuArray(obj.val);
                end

            else

                % Complain and bomb out
                error('incorrect number of input arguments.');

            end

        end

        % RCV sparse matrices are numeric
        function answer=isnumeric(obj) %#ok<MANU>

            % Always yes
            answer=true();

        end

        % RCV sparse matrices are matrices
        function answer=ismatrix(obj) %#ok<MANU>

            % Always yes
            answer=true();

        end

        % RCV sparse matrices are floats
        function answer=isfloat(obj) %#ok<MANU>

            % Always yes
            answer=true();

        end

    end

end

% Consistency enforcement
function grumble(varargin)

% Needs some thought

end

% There are weeks when decades happen.
%
% Vladimir Lenin

