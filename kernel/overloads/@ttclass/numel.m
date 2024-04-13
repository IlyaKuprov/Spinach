% Number of elements in the matrix represented by a tensor 
% train. Syntax:
%
%                         n=numel(tt)
%
% Parameters:
%
%    tt  - tensor train object
%
% Ouputs:
%
%     n  - an integer
%
% Note: for large spin systems, the result may be too large
%       to fit into the 64-bit integer.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/numel.m>

function n=numel(tt)

% Check consistency
grumble(tt);

% Compute the number of elements
n=prod(prod(sizes(tt)));

% Check for infinities
if n>intmax
    error('the number of elements exceeds Matlab''s intmax.');
end

end

% Consistency enforcement
function grumble(tt)
if ~isa(tt,'ttclass')
    error('this function only applies to tensor trains.');
end
end

% If it had been possible to build the tower of Babel without
% ascending it, the work would have been permitted. 
%
% Franz Kafka

