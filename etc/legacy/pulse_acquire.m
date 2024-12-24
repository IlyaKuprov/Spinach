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
% <https://spindynamics.org/wiki/index.php?title=pulse_acquire.m>

function varargout=pulse_acquire(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use acquire(), hp_acquire() or sp_acquire() instead.');

end

% I have brought myself, by long meditation, to the conviction 
% that a human being with a settled purpose must accomplish it, 
% and that nothing can resist a will which would stake even 
% existence upon its fulfillment. 
%
% Benjamin Disraeli

