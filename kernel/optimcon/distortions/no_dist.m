% A distortion function that applies no distortion and therefore
% has a unit Jacobian. Syntax:
%
%                        [w,J]=no_dist(w)
%
% Parameters:
%
%    w         - waveform in rad/s nutation frequency units,
%                one time slice per column, and rows arran-
%                ged as XYXY... with respect to in-phase and
%                quadrature parts on each control channel
%
% Outputs:
%
%    w         - the same waveform as the input
%
%    J         - a sparse unit matrix with the dimension mat-
%                ching the vectorisation of the input
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=no_dist.m>

function [w,J]=no_dist(w)

% Return a unit Jacobian if asked
if nargout>1, J=speye(numel(w)); end

end

% If I only knew how I could get mathematicians interested in
% transformation groups and the treatment of differential equ-
% ations that arises from them. I am absolutely certain, that,
% at some point in the future, these theories will be recogni-
% zed as fundamental.
% 
% Sophus Lie, in his letter
% to Adolph Mayer

