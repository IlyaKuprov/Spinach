% A trap for legacy function calls. At some point, this
% function was a part of Spinach, and may have been men-
% tioned in various published papers.
%
% If it appears in this legacy directory, this function
% was either superceded by somethig more general and po-
% werful, or renamed into something with a more informa-
% tive name, which is printed in the error message.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=expv.m>

function varargout=expv(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use step() instead.');

end

% Of all that is written, I love only what a man 
% has written with his own blood.
%
% Friedrich Nietzsche

