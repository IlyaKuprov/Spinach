% Adds prefix matrices to a polyadic. Anything the polyadic
% multiplies will subsequently be multiplied by the prefix
% matrices. Syntax:
%
%                       p=prefix(a,p)
%
% Parameters:
%
%      a   -   prefix matrix
%
%      p   -   polyadic object
%
% Outputs:
%
%      p   -   polyadic object
%
% Note: a prefix can be a polyadic itself.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/prefix.m>

function p=prefix(a,p)

% Check consistency
grumble(p);

% Represent scalar prefix as a dimensioned scaled identity
if ~isa(a,'polyadic')&&isscalar(a), a=opium(size(p,1),a); end

% Check the dimensions
if size(a,2)~=size(p,1)
    error('matrix dimension mismatch.');
end

% Update prefix array
p.prefix=[{a} p.prefix];

end

% Consistency enforcement
function grumble(p)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
end

% It is difficult to get a man to understand 
% something when his salary depends on his not
% understanding it.
%
% Upton Sinclair


