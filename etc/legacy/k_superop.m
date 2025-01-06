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
% <https://spindynamics.org/wiki/index.php?title=k_superop.m>

function varargout=k_superop(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use kinetics() instead.');

end

% "Why 100? If I were wrong, one would have been enough."
%
% Albert Eistein's response to Nazis lining up
% 100 Aryan scientists to denounce his theories

