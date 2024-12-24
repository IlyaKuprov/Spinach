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
% <https://spindynamics.org/wiki/index.php?title=sim_pulse.m>

function varargout=sim_pulse(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use shaped_pulse_xy() instead.');

end

% Among other evils which being unarmed brings you, 
% it causes you to be despised.
%
% Nicolo Machiavelli

