% Probability density estimation for a three-dimensional Cartesian 
% point cloud on a user-specified regular grid. Syntax:
%
%         density=xyz2pd(coords,x_range,y_range,z_range,...
%                               x_npts, y_npts, z_npts)
%
% Parameters:
%
%    coords  - an N-by-3 array of Cartesian coordinates
%
%    x_range - a two-element vector [xmin xmax] specifying the
%              Cartesian grid extent along the x axis
%
%    y_range - a two-element vector [ymin ymax] specifying the
%              Cartesian grid extent along the y axis
%
%    z_range - a two-element vector [zmin zmax] specifying the
%              Cartesian grid extent along the z axis
%
%    x_npts  - the number of grid points along the x axis
%
%    y_npts  - the number of grid points along the y axis
%
%    z_npts  - the number of grid points along the z axis
%
% Outputs:
%
%    density - a three-dimensional array containing the num-
%              ber of points falling into each grid cell
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=xyz2pd.m>

function density=xyz2pd(coords,x_range,y_range,z_range,...
                               x_npts, y_npts, z_npts)

% Check consistency
grumble(coords,x_range,y_range,z_range,x_npts,y_npts,z_npts);

% Compute grid cell edges
x_edges=linspace(x_range(1),x_range(2),x_npts+1);
y_edges=linspace(y_range(1),y_range(2),y_npts+1);
z_edges=linspace(z_range(1),z_range(2),z_npts+1);

% Assign each point to its grid cell
cell_idx_x=discretize(coords(:,1),x_edges);
cell_idx_y=discretize(coords(:,2),y_edges);
cell_idx_z=discretize(coords(:,3),z_edges);

% Discard points outside the grid
valid_mask=(~isnan(cell_idx_x))&...
           (~isnan(cell_idx_y))&...
           (~isnan(cell_idx_z));
cell_idx_x=cell_idx_x(valid_mask); 
cell_idx_y=cell_idx_y(valid_mask); 
cell_idx_z=cell_idx_z(valid_mask);

% Assign a linear index to each cell
linear_index=sub2ind([x_npts     y_npts     z_npts],...
                      cell_idx_x,cell_idx_y,cell_idx_z);

% Count the number of points falling into each cell
density=accumarray(linear_index,1,[x_npts*y_npts*z_npts 1]);

% Reshape to create a three-dimensional density
density=reshape(density,[x_npts y_npts z_npts]);

end

% Consistency enforcement
function grumble(coords,x_range,y_range,z_range,...
                        x_npts, y_npts, z_npts)
if (~isnumeric(coords))||(size(coords,2)~=3)||...
   (~isreal(coords))||any(~isfinite(coords),'all')
    error('coords must be a real finite matrix with three columns.');
end
if isempty(coords)
    error('coords must contain at least one point.');
end
if (~isnumeric(x_range))||(numel(x_range)~=2)||...
   (~isreal(x_range))||any(~isfinite(x_range))||(x_range(1)>=x_range(2))
    error('x_range must be a two-element real vector with x_range(1)<x_range(2).');
end
if (~isnumeric(y_range))||(numel(y_range)~=2)||...
   (~isreal(y_range))||any(~isfinite(y_range))||(y_range(1)>=y_range(2))
    error('y_range must be a two-element real vector with y_range(1)<y_range(2).');
end
if (~isnumeric(z_range))||(numel(z_range)~=2)||...
   (~isreal(z_range))||any(~isfinite(z_range))||(z_range(1)>=z_range(2))
    error('z_range must be a two-element real vector with z_range(1)<z_range(2).');
end
if (~isnumeric(x_npts))||(~isscalar(x_npts))||...
   (~isreal(x_npts))||(x_npts<2)||(mod(x_npts,1)~=0)
    error('x_npts must be a real integer greater than 1.');
end
if (~isnumeric(y_npts))||(~isscalar(y_npts))||...
   (~isreal(y_npts))||(y_npts<2)||(mod(y_npts,1)~=0)
    error('y_npts must be a real integer greater than 1.');
end
if (~isnumeric(z_npts))||(~isscalar(z_npts))||...
   (~isreal(z_npts))||(z_npts<2)||(mod(z_npts,1)~=0)
    error('z_npts must be a real integer greater than 1.');
end
end

% It's not DNS.
% There's no way it's DNS.
% It was DNS.

