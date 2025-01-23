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
% <https://spindynamics.org/wiki/index.php?title=h_superop.m>

function varargout=h_superop(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use hamiltonian() instead.');

end

% History has shown us that it's not religion that's
% the problem, but any system of thought that insists
% that one group of people are inviolably in the right,
% whereas the others are in the wrong and must somehow
% be punished.
%
% Rod Liddle

