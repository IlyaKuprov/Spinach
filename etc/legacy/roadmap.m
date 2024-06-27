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
% <https://spindynamics.org/wiki/index.php?title=roadmap.m>

function varargout=roadmap(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use powder() instead.');

end

% The fact that we live at the bottom of a deep gravity well, on the
% surface of a gas covered planet going around a nuclear fireball 90
% million miles away and think this to be normal is obviously some
% indication of how skewed our perspective tends to be.
%
% Douglas Adams

