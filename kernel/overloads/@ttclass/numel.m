% Number of elements in the matrix represented by a tensor 
% train. Syntax:
%
%                         n=numel(tt)
%
% Parameters:
%
%    tt  - tensor train object
%
% Outputs:
%
%     n  - an integer
%
% Note: for large spin systems, the result may be too large
%       to be represented exactly as a double.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/numel.m>

function n=numel(tt)

% Check consistency
grumble(tt);

% Compute the number of elements exactly
n=prod(int64(sizes(tt)),'all','native');

% Check for overflow
if n>flintmax
    error('the number of elements exceeds Matlab''s flintmax.');
end

% Return a double
n=double(n);

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

