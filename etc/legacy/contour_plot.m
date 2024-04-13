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
% <https://spindynamics.org/wiki/index.php?title=contour_plot.m>

function varargout=contour_plot(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use plot_2d() instead.');

end

% According to a trade legend, Uhlenbeck and Goudsmit (students of
% Ehrenfest when they stumbled upon the concept of spin) presented
% it to Ehrenfest and said, in effect "here's our theory, but don't
% publish it - it can't be right". He submitted it anyway with the
% justification that they were "young enough to be able to afford
% a stupidity".

