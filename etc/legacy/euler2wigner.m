% A trap for legacy function calls. At some point, this
% function was a part of Spinach, and may have been men-
% tioned in various published papers.
%
% If it appears in this legacy directory, this function
% was either superceded by somethig more general and po-
% werful, or renamed into something with a more informa-
% tive name, which is printed in the error message.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=euler2wigner.m>

function varargout=euler2wigner(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use wigner() instead.');

end

% Always code as if the guy who ends up maintaining 
% your code will be a violent psychopath who knows 
% where you live.
%
% Martin Golding

