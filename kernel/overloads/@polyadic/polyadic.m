% Creates an object of a polyadic class. Syntax:
%
%                          p=polyadic(cores)
%
% A polyadic is a matrix formed by a Kronecker product, with those
% products stored unopened. For example,
%
%                        cores={{A,B,C},{D,E}}
%
% corresponds to A(x)B(x)C + D(x)E matrix. Any multiplicative acti-
% on by this matrix may be computed without opening the Kronecker
% products. This can save orders of magnitude in CPU time. 
%
% Parameters:
%
%         cores -    a cell array of cell arrays of matrices
%                    whose Konecker products make up the mat-
%                    rix of interest.
%
% Outputs:
%
%         p     -    a polyadic representation of a matrix 
%                    that behaves in many respects like the
%                    matrix it represents.
%
% Note: nested polyadics are permitted - the input matrices may be
%       polyadics themselves.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic.m>

classdef (InferiorClasses={?gpuArray,?opium}) polyadic
    
    % Default properties
    properties
        cores={{[]}};
        prefix={};
        suffix={};
    end
    
    % Method description
    methods
        
        % Constructor function
        function p=polyadic(cores)
            
            % Check consistency
            grumble(cores);
            
            % Store the cores
            p.cores=cores;
            
            % Validate the result
            validate(p);
            
        end
        
        % Usually not unit matrix
        function n=iseye(pdc) %#ok<MANU>
            
            n=false(); % Always
            
        end
        
        % Number of elements
        function answer=numel(p)

            % Multiply the dimensions
            answer=prod(size(p)); %#ok<PSIZE>
            
        end
        
        % Polyadics are numeric
        function answer=isnumeric(p) %#ok<MANU>
            
            answer=true(); % Always
            
        end
        
        % Polyadics are matrices
        function answer=ismatrix(p) %#ok<MANU>
            
            answer=true(); % Always
            
        end
        
        % Polyadics are floats
        function answer=isfloat(p) %#ok<MANU>
            
            answer=true(); % Always
            
        end
        
    end
    
end

% Consistency enforcement
function grumble(cores)
if ~iscell(cores)
    error('cores must be a cell array.');
end
for n=1:numel(cores)
    if ~iscell(cores{n})
        error('elements of cores must also be cell arrays.');
    end
end
end

% Sad news arrived from the University of Nottingham as it announced that
% a garden snail whose love life captured the imagination of a nation had
% died. Jeremy -- a "one-in-a-million mutant snail" thanks to his rare le-
% ft-coiling shell -- shot to fame when his handler, Nottingham scientist
% Angus Davison, appeared on BBC Radio 4's Today programme to appeal for a 
% mate for his lonely "leftie" snail. The snail "shellebrity" quickly gai-
% ned cult status: one fan penned a tragic love ballad about Jeremy's pli-
% ght and another had a tattoo of the leftie-coil in his honour. Jeremy's
% burgeoning celebrity, however, appeared to have little effect on poten-
% tial mates. Two potential leftie beaux couriered from Ipswich and Major-
% ca preferred each other and produced more than 300 tiny snail babies be-
% tween them, although Jezza was later to mate with the Spanish import af-
% ter his love rival returned to Suffolk. Shell-shocked fans of Jeremy can
% pay their final respects at the Natural History Museum, where his shell
% will be displayed.
%
% Times Higher Education Magazine, 26 Oct 2017

