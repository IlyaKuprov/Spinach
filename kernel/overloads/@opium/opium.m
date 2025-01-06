% Object Pretending It is a Unit Matrix (OPIUM). Syntax:
%
%                    M=opium(dim,coeff) 
%
% Parameters:
%
%   dim   - dimension of the unit matrix
%
%   coeff - coefficient in front of the unit matrix
%
% Outputs:
%
%   M     - an OPIUM representing the specified matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=opium/opium.m>

classdef (InferiorClasses={?gpuArray}) opium
    
    % Default properties
    properties
        dim=1;
        coeff=1;
    end
    
    % Method description
    methods
        
        % Constructor function
        function M=opium(dim,coeff)
            
            % Check consistency
            grumble(dim,coeff);
            
            % Store the parameters
            M.dim=dim; M.coeff=coeff;
            
        end
        
        % Number of non-zeroes
        function n=nnz(op) %#ok<MANU>
            
            n=1; % Always
            
        end
        
        % Number of elements
        function n=numel(op) %#ok<MANU>
            
            n=1; % Always
            
        end
        
        % Numeric property
        function n=isnumeric(op) %#ok<MANU>
            
            n=true(); % Always
            
        end
        
        % Numeric property
        function n=ismatrix(op) %#ok<MANU>
            
            n=true(); % Always
            
        end
        
        % Is actually unit
        function n=iseye(op)
            
            if op.coeff==1
                n=true();
            else
                n=false();
            end
            
        end
        
        % Conversion to sparse
        function op=sparse(op)
            
            % Return sparse matrix
            op=op.coeff*speye(op.dim);
            
        end
        
        % Conjugation
        function op=conj(op)
            
            % Conjugate the coefficient
            op.coeff=conj(op.coeff);
            
        end
        
        % Conjugate-transpose
        function op=ctranspose(op)
            
            % Conjugate the coefficient
            op.coeff=conj(op.coeff);
            
        end
        
        % GPU upload
        function op=gpuArray(op)

            % Upload coefficient
            op.coeff=gpuArray(op.coeff);

        end
        
        % GPU gather
        function op=gather(op)

            % Gather coefficient
            op.coeff=gather(op.coeff);

        end
        
    end
    
end

% Consistency enforcement
function grumble(dim,coeff)
if (~isnumeric(dim))||(~isreal(dim))||(~isscalar(dim))||...
   (mod(dim,1)~=0)||(dim<1)
    error('dim must be a positive integer scalar.');
end
if (~isnumeric(coeff))||(~isscalar(coeff))
    error('coeff must be a scalar.');
end
end

% If you have a programming question to ask in a forum, ask it
% from a girl's account. You will receive hundreds of replies
% with detailed descriptions of the nature of the problem, and
% of the possible solutions. A few particularly silly wankers
% would even write the code for you.
%
% Russian internet wisdom

