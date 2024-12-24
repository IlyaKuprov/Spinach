% This function packages everything it receives into a cell array and 
% returns it back. This is useful for pulling information back from
% various Spinach wrappers - call this as a pulse sequence. Syntax:
%
%                     answer=impound(varargin)
%
% Parameters:
%
%   varargin   - any number of parameters of any type
%
% Outputs:
%
%   answer     - all input parameters as a cell array
% 
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=impound.m>

function answer=impound(varargin)

% Return what was received
answer=varargin;

end

% He who dares not offend cannot be honest.
%
% Thomas Paine

% #NGRUM