% Switches X axis tick labels to engineering notation by setting
% the numeric-ruler exponent to a multiple of three. Syntax:
%
%                         kxtickfix()
%
% Parameters:
%
%    none
%
% Outputs:
%
%    updates the current axis system
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kxtickfix.m>

function kxtickfix()

% Pull out the current X axis ruler
ax=gca;
ruler=ax.XAxis;

% Check consistency
grumble(ruler);

% Extract finite non-zero tick values
ticks=ruler.TickValues;
ticks=ticks(isfinite(ticks)&(ticks~=0));

% Use finite non-zero limits when tick values are unavailable
if isempty(ticks)
    ticks=ruler.Limits;
    ticks=ticks(isfinite(ticks)&(ticks~=0));
end

% Use unscaled labels for zero axes
if isempty(ticks)
    ruler.Exponent=0;
else
    ruler.Exponent=3*floor(log10(max(abs(ticks)))/3);
end

end

% Consistency enforcement
function grumble(ruler)
if ~isa(ruler,'matlab.graphics.axis.decorator.NumericRuler')
    error('X axis must be numeric.');
end
if ~strcmp(ruler.Scale,'linear')
    error('X axis scale must be linear.');
end
end

