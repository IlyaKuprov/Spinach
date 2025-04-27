% Voronoi tessellation of a 2D COMSOL mesh. Syntax:
%
%                 mesh=mesh_vorn(mesh)
%
% Parameters:
%
%    mesh - Spinach mesh object
%
% Outputs:
%
%    mesh - updated mesh object
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mesh_vorn.m>

function mesh=mesh_vorn(mesh)

% Check consistency
grumble(mesh);

% Run Voronoi tessellation of the mesh
[V,C]=voronoin([mesh.x ...
                mesh.y]);
mesh.vor.vertices=V; 
mesh.vor.cells=C;

% Keep only active cells
mesh.vor.cells=mesh.vor.cells(mesh.idx.active);
mesh.vor.ncells=numel(mesh.vor.cells);

% Voronoi cell area calculation
vor_cell_areas=zeros(mesh.vor.ncells,1);
for n=1:mesh.vor.ncells
    vor_xcoords=mesh.vor.vertices(mesh.vor.cells{n},1);
    vor_ycoords=mesh.vor.vertices(mesh.vor.cells{n},2);
    vor_cell_areas(n)=polyarea(vor_xcoords,vor_ycoords);  
end 

% Add weights to mesh structure 
mesh.vor.weights=vor_cell_areas;

% Find the maximum number of vertices making up the cell
mesh.vor.max_cell_size=max(cellfun(@numel,mesh.vor.cells));

end

% Consistency enforcement
function grumble(mesh)
if ~isfield(mesh,'idx')
    error('vertex index is missing from mesh structure.');
end
end

% A man's rights rest in three boxes: the ballot box,
% the jury box, and the cartridge box.
%
% Frederick Douglass

