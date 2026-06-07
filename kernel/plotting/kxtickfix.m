% Switches X axis tick labels to engineering notation by setting
% the numeric-ruler exponent to a multiple of three and installs
% a pan/zoom auto-updater. Syntax:
%
%                         kxtickfix()
%
% Parameters:
%
%    none
%
% Outputs:
%
%    updates the current axis system and installs an auto-updater
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

% Preserve any existing limit-change callback
app_key='SpinachKXTickFix';
if isappdata(ax,app_key)
    state=getappdata(ax,app_key);
    if isfield(state,'callback')&&isequal(ruler.LimitsChangedFcn,state.callback)
        old_fcn=state.old_fcn;
    else
        old_fcn=ruler.LimitsChangedFcn;
    end
else
    old_fcn=ruler.LimitsChangedFcn;
end

% Install the automatic exponent updater
callback=@(source,event)local_callback(source,event,ax,app_key);
state.callback=callback; state.old_fcn=old_fcn;
setappdata(ax,app_key,state);
ruler.LimitsChangedFcn=callback;

% Set the current exponent
local_update(ruler);

end

% Local update callback
function local_callback(ruler,event,ax,app_key)

% Ignore deleted graphics objects
if ~ishandle(ruler)||~ishandle(ax)||~isappdata(ax,app_key), return, end

% Update the engineering exponent
local_update(ruler);

% Run the pre-existing callback if there was one
state=getappdata(ax,app_key);
if ~isempty(state.old_fcn)
    if isa(state.old_fcn,'function_handle')
        state.old_fcn(ruler,event);
    elseif iscell(state.old_fcn)
        feval(state.old_fcn{1},ruler,event,state.old_fcn{2:end});
    elseif ischar(state.old_fcn)||isstring(state.old_fcn)
        evalin('base',char(state.old_fcn));
    end
end

end

% Local exponent updater
function local_update(ruler)

% Ignore non-linear axis state after installation
if ~strcmp(ruler.Scale,'linear'), return, end

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

