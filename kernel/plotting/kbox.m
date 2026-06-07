% Creates a tickless boxed frame around the current axes 
% using ordinary line objects that live in the same data
% space as the plot. This is needed for spectrograms be-
% cause Matlab has dumb plotting defaults. Syntax:
%
%                       kbox()
%
% Outputs:
%
%    creates or updates a tickless axis box 
%    in the current axes
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kbox.m>

function kbox() % #NGRUM

% Get the current axis object
ax=gca();

% Remove any previous Spinach box overlay
if isappdata(ax,'SpinachKBox')
    state=getappdata(ax,'SpinachKBox'); rmappdata(ax,'SpinachKBox');
    if isfield(state,'listeners')
        cellfun(@delete,state.listeners(cellfun(@isvalid,state.listeners)));
    end
    if isfield(state,'box_line')&&ishandle(state.box_line)
        delete(state.box_line);
    end
    if isfield(state,'box_axes')&&ishandle(state.box_axes)
        delete(state.box_axes);
    end
end

% Delete orphaned overlay axes from the rejected implementation
delete(findall(get(ax,'Parent'),'Type','axes','Tag','SpinachKBox'));

% Create the box line in the plot axes
set(ax,'Box','off','Layer','top');
state.box_line=line(ax,nan,nan,'Color',get(ax,'XColor'),...
                    'LineWidth',get(ax,'LineWidth'),'Clipping','off',...
                    'HandleVisibility','off','HitTest','off',...
                    'PickableParts','none','Tag','SpinachKBox');

% Exclude the box line from autoscaling
lim_props={'XLimInclude','YLimInclude','ZLimInclude'};
for n=1:numel(lim_props)
    if isprop(state.box_line,lim_props{n})
        set(state.box_line,lim_props{n},'off');
    end
end

% Update only on data-space changes, not on camera motion
props={'XLim','YLim','ZLim','LineWidth','XColor','Children'};
state.listeners={};
for n=1:numel(props)
    prop=findprop(ax,props{n});
    if ~isempty(prop)&&prop.SetObservable
        state.listeners{end+1}=addlistener(ax,props{n},'PostSet',...
                                          @(~,~)local_update(ax));
    end
end
state.listeners{end+1}=addlistener(ax,'MarkedClean',@(~,~)local_update(ax));
state.busy=false; setappdata(ax,'SpinachKBox',state);

% Draw the box
local_update(ax);

end

% Local update helper
function local_update(ax)

% Ignore inactive or recursive updates
if ~ishandle(ax)||~isappdata(ax,'SpinachKBox'), return, end
state=getappdata(ax,'SpinachKBox');
if state.busy, return, end
state.busy=true; setappdata(ax,'SpinachKBox',state);

% Build the box line strip
xl=get(ax,'XLim'); yl=get(ax,'YLim'); zl=get(ax,'ZLim');
signature={numel(allchild(ax)),xl,yl,zl,get(ax,'LineWidth'),get(ax,'XColor')};
if isfield(state,'signature')&&isequaln(signature,state.signature)
    state.busy=false; setappdata(ax,'SpinachKBox',state); return
end
state.signature=signature;
is3d=(abs(get(ax,'View')*[0;1]-90)>sqrt(eps))||...
     any(arrayfun(@(h)isprop(h,'ZData')&&~isempty(get(h,'ZData')),allchild(ax)));
x=[]; y=[]; z=[];
if is3d
    for q=yl, for r=zl, [x,y,z]=local_seg(x,y,z,[xl(1) q r],[xl(2) q r]); end, end
    for q=xl, for r=zl, [x,y,z]=local_seg(x,y,z,[q yl(1) r],[q yl(2) r]); end, end
    for q=xl, for r=yl, [x,y,z]=local_seg(x,y,z,[q r zl(1)],[q r zl(2)]); end, end
else
    [x,y,z]=local_seg(x,y,z,[xl(1) yl(1) 0],[xl(2) yl(1) 0]);
    [x,y,z]=local_seg(x,y,z,[xl(2) yl(1) 0],[xl(2) yl(2) 0]);
    [x,y,z]=local_seg(x,y,z,[xl(2) yl(2) 0],[xl(1) yl(2) 0]);
    [x,y,z]=local_seg(x,y,z,[xl(1) yl(2) 0],[xl(1) yl(1) 0]);
end
if is3d
    set(state.box_line,'XData',x,'YData',y,'ZData',z,...
                       'Color',get(ax,'XColor'),'LineWidth',get(ax,'LineWidth'));
else
    set(state.box_line,'XData',x,'YData',y,'ZData',[],...
                       'Color',get(ax,'XColor'),'LineWidth',get(ax,'LineWidth'));
end

% Keep the box above grid and data objects
kids=allchild(ax); box_kid=(kids==state.box_line);
set(ax,'Children',[state.box_line; kids(~box_kid)]);
state.busy=false; setappdata(ax,'SpinachKBox',state);

end

% NaN-separated segment helper
function [x,y,z]=local_seg(x,y,z,a,b)

% Add one segment
x=[x a(1) b(1) nan];
y=[y a(2) b(2) nan];
z=[z a(3) b(3) nan];

end

% "No one realized that the pumps that delivered fuel to the 
% emergency generators were electric."
%
% Angel Feliciano, representative of Verizon
% explaining why Verizon's backup power failed
% during the August 13, 2003 blackout, causing
% disruption to the 911 service.

