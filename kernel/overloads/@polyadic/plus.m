% Polyadic addition operation. Does not perform the actual additi-
% on, but instead stores the operands as a sum of unopened Kronec-
% ker products. Syntax:
%
%                            c=plus(a,b)
%
% Parameters:
%
%   a,b   - polyadic objects
%
% Outputs:
%
%   c     - polyadic object
%
% Note: use this operation sparingly - the additions are simply 
%       buffered, and all subsequent operations will be slower.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/plus.m>

function c=plus(a,b)

% Check consistency
grumble(a,b);

% Run shortcuts
if nnz(a)==0, c=simplify(b); return, end
if nnz(b)==0, c=simplify(a); return, end

% Possible cases
if ~isa(a,'polyadic')
    
    % Matrix + polyadic
    if isempty(b.prefix)&&isempty(b.suffix)
        c=b; c.cores=[c.cores {{a}}];
    else
        c=polyadic({{a},{b}});
    end
        
elseif ~isa(b,'polyadic')
    
    % Polyadic + matrix
    if isempty(a.prefix)&&isempty(a.suffix)
        c=a; c.cores=[c.cores {{b}}];
    else
        c=polyadic({{a},{b}});
    end
    
else
    
    % Polyadic + polyadic 
    if isempty(a.prefix)&&isempty(a.suffix)&&...
       isempty(b.prefix)&&isempty(b.suffix)
       c=polyadic([a.cores b.cores]);
    else
       c=polyadic({{a},{b}});
    end
       
end

% Simplify the result
c=simplify(c);

end

% Consistency enforcement
function grumble(a,b)
[nrows_a,ncols_a]=size(a); 
[nrows_b,ncols_b]=size(b);
if (~isscalar(a))&&(~isscalar(b))
    if (nrows_a~=nrows_b)||(ncols_a~=ncols_b)
        error('operands must represent matrices of the same dimension.');
    end
end
end

% The best revenge is massive success.
%
% Frank Sinatra

