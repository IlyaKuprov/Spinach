% A replacement for the 'grid' command in MATLAB that produces
% grey, rather than black-and-transparent, grid lines suitable
% for publishing. Syntax:
%
%                          kgrid()
%
% Outputs:
%
%    creates or updates a custom grid overlay in the current axes
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kgrid.m>

function kgrid()

% Get the current axis object
ax=gca();

% Remove any previous Spinach grid overlay
if isappdata(ax,'SpinachKGrid')
    state=getappdata(ax,'SpinachKGrid'); rmappdata(ax,'SpinachKGrid');
    if isfield(state,'listeners')
        cellfun(@delete,state.listeners(cellfun(@isvalid,state.listeners)));
    end
    if isfield(state,'minor_line')&&ishandle(state.minor_line)
        delete(state.minor_line);
    end
    if isfield(state,'major_line')&&ishandle(state.major_line)
        delete(state.major_line);
    end
end

% Set native grid styling before replacing it
if ~isappdata(ax,'SpinachKBox'), box(ax,'on'); end
grid(ax,'on'); grid(ax,'minor');
set(ax,'MinorGridAlpha',1,'MinorGridColor',[0.95 0.95 0.95],...
       'MinorGridLineStyle','-','GridAlpha',1,...
       'GridColor',[0.85 0.85 0.85],'GridLineStyle','-',...
       'Layer','top','XGrid','off','YGrid','off','ZGrid','off',...
       'XMinorGrid','off','YMinorGrid','off','ZMinorGrid','off');

% Create two grid line strips in the plot axes
state.minor_line=line(ax,nan,nan,'HandleVisibility','off',...
                      'HitTest','off','PickableParts','none',...
                      'Tag','SpinachKGridMinor');
state.major_line=line(ax,nan,nan,'HandleVisibility','off',...
                      'HitTest','off','PickableParts','none',...
                      'Tag','SpinachKGridMajor');

% Exclude custom grid lines from autoscaling
lim_props={'XLimInclude','YLimInclude','ZLimInclude'};
for n=1:numel(lim_props)
    if isprop(state.minor_line,lim_props{n})
        set([state.minor_line state.major_line],lim_props{n},'off');
    end
end

% Update only on data-space changes, not on camera motion
props={'XLim','YLim','ZLim','XTick','YTick','ZTick',...
       'XMinorTick','YMinorTick','ZMinorTick','GridColor',...
       'MinorGridColor','GridLineStyle','MinorGridLineStyle','LineWidth',...
       'Children'};
state.listeners={};
for n=1:numel(props)
    prop=findprop(ax,props{n});
    if ~isempty(prop)&&prop.SetObservable
        state.listeners{end+1}=addlistener(ax,props{n},'PostSet',...
                                          @(~,~)local_update(ax));
    end
end
state.busy=false;
setappdata(ax,'SpinachKGrid',state);

% Draw the custom grid immediately
local_update(ax);

end


function local_update(ax)

% Ignore inactive or recursive updates
if ~ishandle(ax)||~isappdata(ax,'SpinachKGrid'), return, end
state=getappdata(ax,'SpinachKGrid');
if state.busy, return, end
state.busy=true; setappdata(ax,'SpinachKGrid',state);

% Collect axes geometry
xl=get(ax,'XLim'); yl=get(ax,'YLim'); zl=get(ax,'ZLim');
signature={numel(allchild(ax)),xl,yl,zl,get(ax,'XTick'),get(ax,'YTick'),...
           get(ax,'ZTick'),get(ax,'LineWidth'),get(ax,'GridColor'),...
           get(ax,'MinorGridColor')};
if isfield(state,'signature')&&isequaln(signature,state.signature)
    state.busy=false; setappdata(ax,'SpinachKGrid',state); return
end
state.signature=signature;
is3d=(abs(get(ax,'View')*[0;1]-90)>sqrt(eps))||...
     any(arrayfun(@(h)isprop(h,'ZData')&&~isempty(get(h,'ZData')),allchild(ax)));

% Draw minor and major strips
for k=1:2
    if k==1
        ticks={local_ticks(ax,'x','minor'),local_ticks(ax,'y','minor'),local_ticks(ax,'z','minor')};
        line_obj=state.minor_line; colour=get(ax,'MinorGridColor'); style=get(ax,'MinorGridLineStyle');
    else
        ticks={local_ticks(ax,'x','major'),local_ticks(ax,'y','major'),local_ticks(ax,'z','major')};
        line_obj=state.major_line; colour=get(ax,'GridColor'); style=get(ax,'GridLineStyle');
    end
    x=[]; y=[]; z=[];
    if is3d
        for q=ticks{1}, [x,y,z]=local_seg(x,y,z,[q yl(1) zl(1)],[q yl(2) zl(1)]); end
        for q=ticks{2}, [x,y,z]=local_seg(x,y,z,[xl(1) q zl(1)],[xl(2) q zl(1)]); end
        for q=ticks{1}, [x,y,z]=local_seg(x,y,z,[q yl(1) zl(1)],[q yl(1) zl(2)]); end
        for q=ticks{3}, [x,y,z]=local_seg(x,y,z,[xl(1) yl(1) q],[xl(2) yl(1) q]); end
        for q=ticks{2}, [x,y,z]=local_seg(x,y,z,[xl(1) q zl(1)],[xl(1) q zl(2)]); end
        for q=ticks{3}, [x,y,z]=local_seg(x,y,z,[xl(1) yl(1) q],[xl(1) yl(2) q]); end
    else
        for q=ticks{1}, x=[x q q nan]; y=[y yl nan]; end %#ok<AGROW>
        for q=ticks{2}, x=[x xl nan]; y=[y q q nan]; end %#ok<AGROW>
        z=[];
    end
    set(line_obj,'XData',x,'YData',y,'ZData',z,'Color',colour,...
                 'LineStyle',style,'LineWidth',get(ax,'LineWidth'));
end

% Stack data above major grid, and major grid above minor grid
kids=allchild(ax);
grid_kids=(kids==state.minor_line)|(kids==state.major_line);
set(ax,'Children',[kids(~grid_kids); state.major_line; state.minor_line]);
state.busy=false; setappdata(ax,'SpinachKGrid',state);

end


function ticks=local_ticks(ax,dim,kind)

% Read major or minor tick locations
lim=get(ax,[upper(dim) 'Lim']);
major=get(ax,[upper(dim) 'Tick']);
if strcmp(kind,'major')
    ticks=major;
elseif strcmp(get(ax,[upper(dim) 'MinorTick']),'on')
    ruler=get(ax,[upper(dim) 'Axis']);
    ticks=get(ruler,'MinorTickValues');
    for n=1:numel(major)
        ticks=ticks(abs(ticks-major(n))>1024*eps(max([1 abs(lim) abs(major)])));
    end
else
    ticks=[];
end

% Keep only finite interior ticks
tol=1024*eps(max([1 abs(lim)]));
ticks=double(ticks(:)).';
ticks=ticks(isfinite(ticks)&(ticks>min(lim)+tol)&(ticks<max(lim)-tol));

end


function [x,y,z]=local_seg(x,y,z,a,b)

% Add one NaN-separated segment
x=[x a(1) b(1) nan];
y=[y a(2) b(2) nan];
z=[z a(3) b(3) nan];

end


% Frigid gentlewomen of the jury! I had thought that 
% months, perhaps years, would elapse before I dared
% to reveal myself to Dolores Haze; but by six she 
% was wide awake, and by six fifteen we were techni-
% cally lovers. I am going to tell you something very
% strange: it was she who seduced me.
%
% Vladimir Nabokov, "Lolita"

% #NGRUM
