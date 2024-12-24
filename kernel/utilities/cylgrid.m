% Draws a cylindrical grid with 10% spacing added around 
% the indicated data extent values. Syntax:
%
%                  cylgrid(zmin,zmax,rmax)
%
% Parameters:
%
%    zmin - lower bound on the Z axis
%
%    zmax - upper bound on the Z axis
%
%    rmax - upper bound on the radius
%
% Outputs:
%
%    this function updates the current figure
%
% ilya.kuprov@weizmann.ac.uk
% emile.ottoy@ugent.be
%
% <https://spindynamics.org/wiki/index.php?title=cylgrid.m>

function cylgrid(zmin,zmax,rmax)

% Check consistency
grumble(zmin,zmax,rmax);

% Get the extent gaps
rgap=0.1*rmax; zgap=0.1*(zmax-zmin);

% Draw the spokes
for n=0:30:330
    line([0 (rmax+rgap)*cosd(n)],...
         [0 (rmax+rgap)*sind(n)],...
         [(zmin-zgap) (zmin-zgap)],'Color',[0.8 0.8 0.8]);
    text((rmax+2*rgap)*cosd(n),...
         (rmax+2*rgap)*sind(n),...
         (zmin-zgap),num2str(n),...
         'HorizontalAlignment','center','VerticalAlignment','middle');
    line([0 (rmax+rgap)*cosd(n)],...
         [0 (rmax+rgap)*sind(n)],...
         [(zmax+zgap) (zmax+zgap)],'Color',[0.8 0.8 0.8]);
    line([(rmax+rgap)*cosd(n) (rmax+rgap)*cosd(n)],...
         [(rmax+rgap)*sind(n) (rmax+rgap)*sind(n)],...
         [(zmin-zgap) (zmax+zgap)],'Color',[0.8 0.8 0.8]);
end

% Draw the circles
phi_grid=linspace(0,2*pi,360);
for r=linspace(0,rmax+rgap,7)
    line(r*cos(phi_grid),r*sin(phi_grid),zmin*ones(1,360)-zgap,'Color',[0.8 0.8 0.8]);
    line(r*cos(phi_grid),r*sin(phi_grid),zmax*ones(1,360)+zgap,'Color',[0.8 0.8 0.8]);
end

% Set axis extents
axis([(-rmax-2*rgap) (rmax+2*rgap)...
      (-rmax-2*rgap) (rmax+2*rgap)...
      (+zmin-2*zgap) (zmax+2*zgap)]);
  
% Draw vertical tick labels
tick_vals=linspace(zmin-zgap,zmax+zgap,7);
for n=1:numel(tick_vals)
    line([(rmax+rgap) (rmax+rgap)],[-rmax/50 rmax/50],...
         [tick_vals(n) tick_vals(n)],'Color',[0.8 0.8 0.8]);
    text((rmax+rgap),0,tick_vals(n),[' ' num2str(tick_vals(n))],...
         'HorizontalAlignment','left','VerticalAlignment','middle');
end

% Switch on the perspective, kill default axes
set(gca,'Projection','perspective','Box','off',...
    'XTick',[],'YTick',[],'ZTick',[],'visible','off');
set(gcf,'Color',[1 1 1]);

end

% Consistency enforcement
function grumble(zmin,zmax,rmax)
if (~isnumeric(zmin))||(~isreal(zmin))||...
   (numel(zmin)~=1)||(~isfinite(zmin))
    error('zmin must be a real number.');
end
if (~isnumeric(zmax))||(~isreal(zmax))||...
   (numel(zmax)~=1)||(~isfinite(zmax))
    error('zmax must be a real number.');
end
if zmin>=zmax
    error('zmin must be smaller than zmax.');
end
if (~isnumeric(rmax))||(~isreal(rmax))||...
   (numel(rmax)~=1)||(~isfinite(rmax))||(rmax<=0)
    error('rmax must be a positive real number.');
end
end

% Socialism in general has a record of failure so blatant 
% that only an intellectual could ignore or evade it.
%
% Thomas Sowell

