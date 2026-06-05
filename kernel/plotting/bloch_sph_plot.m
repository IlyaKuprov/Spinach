% Bloch sphere plot of a 3D trajectory with linear inter-
% polation between points. Also plots the normalised tip
% trajectory of the instantaneous rotation axis, or writes
% a trajectory movie with a moving axis stick. Syntax:
%
%              bloch_sph_plot(x,y,z,varargin)
%              bloch_sph_plot(x,y,z,varargin,file_name)
%
% Parameters:
%
%   x,y,z    - arrays of equal size containing
%              trajectory or trajectories (rows)
%
%   varargin - further arguments to be supplied 
%              to the plot3() functon of Matlab
%
%   file_name - optional character string ending in
%               .mp4 or .avi; must be the last argument
%
% Output:
%
%   this function produces a figure; if a file name
%   is supplied, it also writes a movie file
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bloch_sph_plot.m>

function bloch_sph_plot(x,y,z,varargin)

% Look for a movie file name
file_name='';
if (~isempty(varargin))&&ischar(varargin{end})&&(numel(varargin{end})>=4)
    file_ext=varargin{end}((end-3):end);
    if strcmpi(file_ext,'.mp4')||strcmpi(file_ext,'.avi')
        file_name=varargin{end};
        varargin=varargin(1:(end-1));
    end
end

% Check consistency
grumble(x,y,z,file_name);

% Kill stupid ass figure defaults in R2025a and later 
set(groot,'defaultFigurePosition',[680 458 560 420]); 
set(groot,'defaultFigureWindowStyle','normal'); 
set(groot,'defaultFigureMenuBar','figure'); 
set(groot,'defaultFigureToolbar','figure'); 

% Reference unit sphere
kfigure(); [X,Y,Z]=sphere; 
surf(X,Y,Z,'FaceAlpha',0.25,'EdgeAlpha',0.25);
colormap bone; box on; axis square; hold on;
set(gca,'Projection','perspective'); kgrid;

% Finite differences along each trajectory
dx_dt=zeros(size(x)); dy_dt=zeros(size(y)); dz_dt=zeros(size(z));
if size(x,2)>1
    dx_dt(:,1)=x(:,2)-x(:,1);
    dy_dt(:,1)=y(:,2)-y(:,1);
    dz_dt(:,1)=z(:,2)-z(:,1);
    dx_dt(:,end)=x(:,end)-x(:,end-1);
    dy_dt(:,end)=y(:,end)-y(:,end-1);
    dz_dt(:,end)=z(:,end)-z(:,end-1);
end
if size(x,2)>2
    dx_dt(:,2:(end-1))=(x(:,3:end)-x(:,1:(end-2)))/2;
    dy_dt(:,2:(end-1))=(y(:,3:end)-y(:,1:(end-2)))/2;
    dz_dt(:,2:(end-1))=(z(:,3:end)-z(:,1:(end-2)))/2;
end

% Minimal instantaneous rotation axis direction
axis_x=y.*dz_dt-z.*dy_dt;
axis_y=z.*dx_dt-x.*dz_dt;
axis_z=x.*dy_dt-y.*dx_dt;

% Normalise the axis curve into the unit sphere
axis_norm=sqrt(axis_x.^2+axis_y.^2+axis_z.^2);
good_axis=(axis_norm>0);
axis_x(good_axis)=axis_x(good_axis)./axis_norm(good_axis);
axis_y(good_axis)=axis_y(good_axis)./axis_norm(good_axis);
axis_z(good_axis)=axis_z(good_axis)./axis_norm(good_axis);
axis_x(~good_axis)=NaN;
axis_y(~good_axis)=NaN;
axis_z(~good_axis)=NaN;

% Static trajectory and instantaneous axis-tip curve
if isempty(file_name)

    % Pack rows into NaN-separated plot streams
    traj_x=[x NaN(size(x,1),1)]';
    traj_y=[y NaN(size(y,1),1)]';
    traj_z=[z NaN(size(z,1),1)]';
    axis_curve_x=[axis_x NaN(size(axis_x,1),1)]';
    axis_curve_y=[axis_y NaN(size(axis_y,1),1)]';
    axis_curve_z=[axis_z NaN(size(axis_z,1),1)]';

    % Plot the trajectory and instantaneous axis curve
    plot3(traj_x(:),traj_y(:),traj_z(:),varargin{:});
    plot3(axis_curve_x(:),axis_curve_y(:),axis_curve_z(:),...
          'r--','LineWidth',1.5,...
          'Tag','InstantaneousAxis');

else

    % Animated trajectory and instantaneous axis stick
    axis vis3d;
    if strcmpi(file_name((end-3):end),'.mp4')
        writerObj=VideoWriter(file_name,'MPEG-4');
    else
        writerObj=VideoWriter(file_name,'Motion JPEG AVI');
    end
    writerObj.Quality=100; open(writerObj);
    traj_x=[x(:,1) NaN(size(x,1),1)]';
    traj_y=[y(:,1) NaN(size(y,1),1)]';
    traj_z=[z(:,1) NaN(size(z,1),1)]';
    axis_stick_x=[zeros(size(axis_x,1),1) axis_x(:,1) NaN(size(axis_x,1),1)]';
    axis_stick_y=[zeros(size(axis_y,1),1) axis_y(:,1) NaN(size(axis_y,1),1)]';
    axis_stick_z=[zeros(size(axis_z,1),1) axis_z(:,1) NaN(size(axis_z,1),1)]';
    traj_obj=plot3(traj_x(:),traj_y(:),traj_z(:),varargin{:});
    axis_obj=plot3(axis_stick_x(:),axis_stick_y(:),axis_stick_z(:),...
                   'r-','LineWidth',1.5,'Tag','InstantaneousAxis');
    for n=1:size(x,2)

        % Update the magnetisation trajectory
        traj_x=[x(:,1:n) NaN(size(x,1),1)]';
        traj_y=[y(:,1:n) NaN(size(y,1),1)]';
        traj_z=[z(:,1:n) NaN(size(z,1),1)]';
        set(traj_obj,'XData',traj_x(:),'YData',traj_y(:),'ZData',traj_z(:));

        % Update the instantaneous rotation-axis stick
        axis_stick_x=[zeros(size(axis_x,1),1) axis_x(:,n) NaN(size(axis_x,1),1)]';
        axis_stick_y=[zeros(size(axis_y,1),1) axis_y(:,n) NaN(size(axis_y,1),1)]';
        axis_stick_z=[zeros(size(axis_z,1),1) axis_z(:,n) NaN(size(axis_z,1),1)]';
        set(axis_obj,'XData',axis_stick_x(:),...
                     'YData',axis_stick_y(:),...
                     'ZData',axis_stick_z(:));

        % Grab the frame
        drawnow; writeVideo(writerObj,getframe(gcf));

    end
    close(writerObj);
end

end

% Consistency enforcement
function grumble(x,y,z,file_name)
if (~isnumeric(x))||(~isreal(x))||(~ismatrix(x))||any(~isfinite(x),'all')
    error('x must be a finite real matrix.');
end
if (~isnumeric(y))||(~isreal(y))||(~ismatrix(y))||any(~isfinite(y),'all')
    error('y must be a finite real matrix.');
end
if (~isnumeric(z))||(~isreal(z))||(~ismatrix(z))||any(~isfinite(z),'all')
    error('z must be a finite real matrix.');
end
if (~isequal(size(x),size(y)))||(~isequal(size(x),size(z)))
    error('x, y, and z must have the same dimensions.');
end
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% The Proton VPN servers for Afghanistan running at 
% over 80% load was not something that we expected
% to see when we set them up last year. 
%
% David Peterson, after the UK started
% requiring proof of age on porn sites
% in August 2025
