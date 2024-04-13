% Creates an object of a tensor train class. A tensor train is a type
% of un-opened Kronecker product that behaves as a matrix or a vector
% of a very large dimension, but takes a reasonable amount of memory
% to store. See https://doi.org/10.1137/090752286 for further informa-
% tion. Syntax:
%
%                 tt=ttclass(coeff,kronterms,tolerance)
%
% Parameters:
%
%        coeff   - coefficient in front of the spin operator,
%                  usually the interaction magnitude
%
%    kronterms   - column cell array of matrices whose Krone-
%                  cker product makes up the spin operator
%
%    tolerance   - maximum deviation in the 2-norm between the
%                  TT representation and the flat matrix repre-
%                  sentation that the TT format is allowed to
%                  introduce
%
% Outputs:
%
%           tt   - tensor train object
%
% Notes:
%
% 1. If multiple columns are supplied in kronterms, multiple coeffi-
%    cients are given in coeff, and multiple tolerances are given in
%    tolerance, the resulting tensor train is assumed to be the sum
%    of the individual tensor trains specified in different columns.
%
% 2. Tensor trains are exotic and capricious structures, do not use
%    them unless you know what you are doing.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass.m>

classdef ttclass
    
    % Default properties
    properties
        coeff=0;      % Physical coefficients      (row vector)
        cores={0};    % Cores of each tensor train (in columns)
        tolerance=0;  % Accuracy tolerances        (row vector)
        debuglevel=0; % Diagnostic output switch
    end
    
    % Method description
    methods
        
        % Constructor function
        function tt=ttclass(coeff,kronterms,tolerance)
            
            % Check the number of inputs
            if nargin==0, return;
            elseif (nargin==1)||(nargin==2)||(nargin>3)
                error('incorrect number of input arguments.');
            end
                        
            % Check consistency
            grumble(coeff,kronterms,tolerance);
            
            % Convert Kronecker products into tensor trains
            tt.cores=kronterms;
            for n=1:size(kronterms,2)
                for d=1:size(kronterms,1)
                    tt.cores{d,n}=reshape(full(kronterms{d,n}),[1 size(kronterms{d,n}) 1]);
                end
            end
            
            % Store the coefficients
            tt.coeff=coeff;
            
            % Store the tolerances
            tt.tolerance=tolerance;
            
            % Do not print diagnostics
            tt.debuglevel=0;
            
        end
        
        % Number of cores
        function answer=ncores(obj)
            answer=size(obj.cores,1);
        end
        
        % Number of trains in a sum
        function answer=ntrains(obj)
            answer=size(obj.cores,2);
        end
        
        % Matlab subsref bug workaround
        function n=numArgumentsFromSubscript(obj,s,ic)
            n=builtin('numArgumentsFromSubscript',obj,s,ic);
        end
        
    end
    
end

% Consistency enforcement
function grumble(coeff,cores,tolerance)
if (~isnumeric(coeff))||(~isrow(coeff))
    error('coeff must be a complex row vector.');
elseif (~iscell(cores))||isempty(cores)||numel(size(cores))~=2
    error('cores must be a two-dimensional cell array with at least one element.');
elseif any(~cellfun(@isnumeric,cores(:)))||...
       any(~cellfun(@ismatrix,cores(:)))
    error('all elements of the cores cell array must be matrices.');
elseif (~isnumeric(tolerance))||(~isreal(tolerance))||...
       (~isrow(tolerance))||any(tolerance<0)
    error('tolerance must be a row vector with real non-negative elements.');
elseif numel(coeff)~=numel(tolerance)
    error('coefficient and tolerance vectors must have the same number of elements.');
elseif size(cores,2)~=numel(coeff)
    error('the number of columns in the operand array must match the number of coefficients.')
end
end

% "Kuprov cocktail": caffeine + modafinil + propranolol, 50 mg each.
% "Double Kuprov": 100 mg variety, reserved for extreme workloads.

