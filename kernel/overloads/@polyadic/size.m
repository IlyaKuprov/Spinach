% Returns the size of the matrix represented by the polyadic. Syntax: 
%
%                            answer=size(p,dim)
%
% Parameters:
%
%    p   - a polyadic object
%
%    dim - dimension whose size is required
%
% Outputs:
%
%    answer - a vector with one or two elements
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/size.m>

function varargout=size(p,dim)

% Check consistency
if nargin==2, grumble(dim); end

% Get row dimension
if ~isempty(p.prefix)
    
    % The leftmost matrix in the prefix
    nrows=size(p.prefix{1},1);
    
else
    
    % The cores of the polyadic
    nrows=prod(cellfun(@(x)size(x,1),p.cores{1}));
    
end

% Get column dimension
if ~isempty(p.suffix)
    
    % The rightmost matrix in the suffix
    ncols=size(p.suffix{end},2);
    
else
    
    % The cores of the polyadic
    ncols=prod(cellfun(@(x)size(x,2),p.cores{1}));
    
end

% Compose the answer
if (nargin==1)&&(nargout<=1)
    varargout{1}=[nrows ncols];
elseif (nargin==1)&&(nargout==2)
    varargout{1}=nrows;
    varargout{2}=ncols;
elseif (nargin==2)&&(dim==1)
    varargout{1}=nrows;
elseif (nargin==2)&&(dim==2)
    varargout{1}=ncols;
else
    error('invalid call syntax.');
end

end

% Consistency enforcement
function grumble(dim)
if (~isscalar(dim))||(~ismember(dim,[1 2]))
    error('for a polyadic object, dim must be 1 or 2');
end
end

% What is true, just and beautiful is not determined by popular 
% vote. The masses everywhere are ignorant, short-sighted, moti-
% vated by envy, and easy to fool. Democratic politicians must 
% appeal to these masses in order to be elected. Whoever is the
% best demagogue will win. Almost by necessity, then, democracy
% will lead to the perversion of truth, justice, and beauty.
%
% Hans-Hermann Hoppe

