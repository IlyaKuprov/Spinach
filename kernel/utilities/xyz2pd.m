% Probability density estimation for a three-dimensional Cartesian point
% cloud on a user-specified regular grid. Syntax:
%
%  [density,x_grid,y_grid,z_grid]=xyz2pd(coords,x_lim,y_lim,z_lim,...
%                                        x_points,y_points,z_points)
%
% Parameters:
%
%    coords     - an N-by-3 array of Cartesian coordinates
%
%    x_lim      - a two-element vector [xmin xmax] specifying the
%                 Cartesian grid extent along the x axis
%
%    y_lim      - a two-element vector [ymin ymax] specifying the
%                 Cartesian grid extent along the y axis
%
%    z_lim      - a two-element vector [zmin zmax] specifying the
%                 Cartesian grid extent along the z axis
%
%    x_points   - the number of grid points along the x axis
%
%    y_points   - the number of grid points along the y axis
%
%    z_points   - the number of grid points along the z axis
%
% Outputs:
%
%    density    - a three-dimensional array containing the probability
%                 density evaluated on the grid specified
%
%    x_grid     - the x axis values corresponding to the first dimension
%                 of the density array
%
%    y_grid     - the y axis values corresponding to the second dimension
%                 of the density array
%
%    z_grid     - the z axis values corresponding to the third dimension
%                 of the density array
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=xyz2pd.m>

function [density,x_grid,y_grid,z_grid]=xyz2pd(coords,x_lim,y_lim,z_lim,...
                                               x_points,y_points,z_points)

% Check consistency
grumble(coords,x_lim,y_lim,z_lim,...
        x_points,y_points,z_points);

% Compute grid edges
x_edges=linspace(x_lim(1),x_lim(2),x_points+1);
y_edges=linspace(y_lim(1),y_lim(2),y_points+1);
z_edges=linspace(z_lim(1),z_lim(2),z_points+1);

% Compute grid centers
x_grid=0.5*(x_edges(1:end-1)+x_edges(2:end));
y_grid=0.5*(y_edges(1:end-1)+y_edges(2:end));
z_grid=0.5*(z_edges(1:end-1)+z_edges(2:end));

% Get the probability density
bin_vol=(x_edges(2)-x_edges(1))*...
        (y_edges(2)-y_edges(1))*...
        (z_edges(2)-z_edges(1));    % Bin volume
bin_x=discretize(coords(:,1),x_edges);
bin_y=discretize(coords(:,2),y_edges);
bin_z=discretize(coords(:,3),z_edges);
valid_mask=~isnan(bin_x)&~isnan(bin_y)&~isnan(bin_z);
linear_index=sub2ind([x_points,y_points,z_points],bin_x(valid_mask),bin_y(valid_mask),bin_z(valid_mask));
density=accumarray(linear_index,1,[x_points*y_points*z_points,1]);
density=reshape(density,[x_points,y_points,z_points]);
if any(density,'all')
    density=density/(sum(density,'all')*bin_vol);
end

end

% Consistency enforcement
function grumble(coords,x_lim,y_lim,z_lim,...
                 x_points,y_points,z_points)
if (~isnumeric(coords))||(size(coords,2)~=3)||...
   (~isreal(coords))||any(~isfinite(coords(:)))
    error('coords must be a real finite matrix with three columns.');
end
if isempty(coords)
    error('coords must contain at least one point.');
end
if (~isnumeric(x_lim))||(numel(x_lim)~=2)||...
   (~isreal(x_lim))||any(~isfinite(x_lim))||(x_lim(1)>=x_lim(2))
    error('x_lim must be a two-element real vector with x_lim(1)<x_lim(2).');
end
if (~isnumeric(y_lim))||(numel(y_lim)~=2)||...
   (~isreal(y_lim))||any(~isfinite(y_lim))||(y_lim(1)>=y_lim(2))
    error('y_lim must be a two-element real vector with y_lim(1)<y_lim(2).');
end
if (~isnumeric(z_lim))||(numel(z_lim)~=2)||...
   (~isreal(z_lim))||any(~isfinite(z_lim))||(z_lim(1)>=z_lim(2))
    error('z_lim must be a two-element real vector with z_lim(1)<z_lim(2).');
end
if (~isnumeric(x_points))||(~isscalar(x_points))||(~isreal(x_points))||...
   (x_points<2)||(mod(x_points,1)~=0)
    error('x_points must be an integer greater than one.');
end
if (~isnumeric(y_points))||(~isscalar(y_points))||(~isreal(y_points))||...
   (y_points<2)||(mod(y_points,1)~=0)
    error('y_points must be an integer greater than one.');
end
if (~isnumeric(z_points))||(~isscalar(z_points))||(~isreal(z_points))||...
   (z_points<2)||(mod(z_points,1)~=0)
    error('z_points must be an integer greater than one.');
end
end

% It's not DNS.
% There's no way it's DNS.
% It was DNS.

