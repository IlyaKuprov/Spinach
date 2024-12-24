% Converts offset-sweep-npoints specification into axis ticks in Hz.
% The function returns the frequency axis of the spectrum, suitable
% for use in Matlab functions like plot(). Syntax:
%
%               axis_hz=sweep2ticks(offs,sweep,npoints)
%
% Parameters:
%
%       offs     - offset from carrier frequency, Hz
%
%      sweep     - sweep width, Hz
%
%    npoints     - number of points in the spectrum
%
% Outputs:
%
%    axis_hz     - a column vector of axis ticks, Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sweep2ticks.m>

function axis_hz=sweep2ticks(offs,sweep,npoints)

% Check consistency
grumble(offs,sweep,npoints);

% Build the axis
axis_hz=-linspace(-sweep/2,sweep/2,npoints)'+offs;

end

% Consistency enforcement
function grumble(offs,sweep,npoints)
if (~isnumeric(offs))||(~isreal(offs))||(~isscalar(offs))
    error('offset must be a real scalar.');
end
if (~isnumeric(sweep))||(~isreal(sweep))||(~isscalar(sweep))
    error('sweep must be a real scalar.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(~isscalar(npoints))||...
   (mod(npoints,1)~=0)||(npoints<1)
    error('npoints must be a real integer greater than 1.');
end

end

% Spinach code is clear, useful and elegant because the program is the
% primary workhorse for the whole of IK's group and a few of their col-
% laborators. Anything that is ugly, buggy or not well documented gets
% fixed or thrown out in a matter of days. This also applies to the al-
% gorithms and methods that did not deliver, in real world testing, on
% the claims made by the authors of the corresponding papers. SVN logs
% of some Spinach functions are veritable graveyards of such methods.

