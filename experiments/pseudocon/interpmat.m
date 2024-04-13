% Returns a matrix that acts on a stretched pseudocontact shift density cube
% and projects out the values of the PCS at the Cartesian coordinates given.
% Tricubic interpolation is used. Syntax:
%
%                     P=interpmat(cube_dims,ranges,xyz)
%
% Parameters: 
%
%     cube_dims - pseudocontact shift cube grid sizes, a vector of 
%                 three integers ordered as [X Y Z]
%
%     ranges    - cartesian axis extents for the pseudocontact shift 
%                 cube as [xmin xmax ymin ymax zmin zmax] in Angstroms.
%
%     xyz       - nuclear coordinates as [x y z] with multiple rows) at
%                 which PCS is to be evaluated, in Angstroms.
%
% Output:
% 
%     P         - matrix projecting out PCS values at the specified 
%                 nuclear positions from the stretched PCS cube.
%
% Note: this function is a part of the PCS inverse problem solver module; it
%       should not normally be called directly by the user.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=interpmat.m>

function P=interpmat(cube_dims,ranges,xyz)

% Check consistency
grumble(cube_dims,ranges,xyz);

% Extract the bounds
xmin=ranges(1); xmax=ranges(2);
ymin=ranges(3); ymax=ranges(4);
zmin=ranges(5); zmax=ranges(6);

% Compute grids
x_grid=linspace(xmin,xmax,cube_dims(1));
y_grid=linspace(ymin,ymax,cube_dims(2));
z_grid=linspace(zmin,zmax,cube_dims(3));

% Compute grid intervals
x_grid_interval=(xmax-xmin)/(cube_dims(1)-1);
y_grid_interval=(ymax-ymin)/(cube_dims(2)-1);
z_grid_interval=(zmax-zmin)/(cube_dims(3)-1);

% Preallocate return arrays
rows=cell(size(xyz,1),1);
cols=cell(size(xyz,1),1);
vals=cell(size(xyz,1),1);

% Loop over points
parfor n=1:size(xyz,1) %#ok<*PFBNS>

    % Move into fractional coordinates
    x=(xyz(n,1)-xmin)/x_grid_interval;
    y=(xyz(n,2)-ymin)/y_grid_interval;
    z=(xyz(n,3)-zmin)/z_grid_interval;
    
    % Decide stencil points
    x_stencil=[min([max([1 (ceil(x)-1)]) (cube_dims(1)-3)])...
               min([max([2 (ceil(x)-0)]) (cube_dims(1)-2)])... 
               min([max([3 (ceil(x)+1)]) (cube_dims(1)-1)])...
               min([max([4 (ceil(x)+2)]) (cube_dims(1)-0)])];
    y_stencil=[min([max([1 (ceil(y)-1)]) (cube_dims(2)-3)])...
               min([max([2 (ceil(y)-0)]) (cube_dims(2)-2)])... 
               min([max([3 (ceil(y)+1)]) (cube_dims(2)-1)])...
               min([max([4 (ceil(y)+2)]) (cube_dims(2)-0)])];
    z_stencil=[min([max([1 (ceil(z)-1)]) (cube_dims(3)-3)])...
               min([max([2 (ceil(z)-0)]) (cube_dims(3)-2)])... 
               min([max([3 (ceil(z)+1)]) (cube_dims(3)-1)])...
               min([max([4 (ceil(z)+2)]) (cube_dims(3)-0)])];
    
    % Extract subgrid values
    x_subgrid=x_grid(x_stencil); y_subgrid=y_grid(y_stencil); z_subgrid=z_grid(z_stencil); 
           
    % Compute interpolation vectors
    x_intvec=spalloc(1,cube_dims(1),4); x_intvec(x_stencil)=fdweights(xyz(n,1),x_subgrid,0); %#ok<SPRIX>
    y_intvec=spalloc(1,cube_dims(2),4); y_intvec(y_stencil)=fdweights(xyz(n,2),y_subgrid,0); %#ok<SPRIX>
    z_intvec=spalloc(1,cube_dims(3),4); z_intvec(z_stencil)=fdweights(xyz(n,3),z_subgrid,0); %#ok<SPRIX>
   
    % Store non-zeroes
    [~,cols{n},vals{n}]=find(kron(z_intvec,...
                             kron(y_intvec,x_intvec)));
    rows{n}=n*ones(size(cols{n}));
   
end

% Merge cells
rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals); 

% Form the matrix
P=sparse(rows,cols,vals,size(xyz,1),prod(cube_dims));

end

% Consistency enforcement
function grumble(cube_dims,ranges,xyz)
if (~isnumeric(cube_dims))||(~isreal(cube_dims))||...
   (numel(cube_dims)~=3)||any(mod(cube_dims,1)~=0)
    error('cube_dims must contain three positive integer elements.');
end
if any(cube_dims<4), error('too few points in the cube.'); end
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)
    error('xyz must be an Nx3 array of atomic coordinates.');
end
if (~isnumeric(ranges))||(~isreal(ranges))||(numel(ranges)~=6)
    error('ranges must be a real vector with six elements.');
end
if (ranges(1)>=ranges(2))||(ranges(3)>=ranges(4))||(ranges(5)>=ranges(6))
    error('ranges array should have xmin<xmax, ymin<ymax and zmin<zmax.');
end
if any(xyz(:,1)<ranges(1))||any(xyz(:,1)>ranges(2))||...
   any(xyz(:,2)<ranges(3))||any(xyz(:,2)>ranges(4))||...
   any(xyz(:,3)<ranges(5))||any(xyz(:,3)>ranges(6))
    error('extrapolation is not supported.');
end
end

% Pick up a camera. Shoot something. No matter how small, no matter how
% cheesy, no matter whether your friends and your sister star in it. Put
% your name on it as director. Now you're a director. Everything after
% that - you're just negotiating your budget and your fee.
%
% James Cameron

