% Returns true for polyadics that represent a matrix with a 
% zero dimension. Syntax:
%
%                      answer=isempty(p)
%
% Parameters:
%
%       p  - a polyadic object
%
% Outputs:
%
%   answer - a logical value
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/isempty.m>

function answer=isempty(p)

% Check the size
if any(size(p)==0) 
    answer=true();
else
    answer=false();
end

end

% Q: Why do sumo wrestlers shave their legs and armpits?
% A: To make sure people can tell them apart from feminists.
%
% A "festive season" cracker

