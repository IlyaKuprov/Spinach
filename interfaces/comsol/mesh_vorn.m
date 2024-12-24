% Voronoi tessellation of a 2D COMSOL mesh. Syntax:
%
%        spin_system=mesh_vorn(spin_system)
%
% Parameters:
%
%    spin_system - Spinach data structure with a .mesh
%                  subfield present
%
% Outputs:
%
%    spin_system - Spinach data structure with a .vor
%                  subfield added to the mesh structure
%
% ilya.kuprov@weizmann.ac.uk
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mesh_vorn.m>

function spin_system=mesh_vorn(spin_system)

% Check consistency
grumble(spin_system);

% Run Voronoi tessellation of the mesh
[V,C]=voronoin([spin_system.mesh.x ...
                spin_system.mesh.y]);
spin_system.mesh.vor.vertices=V; 
spin_system.mesh.vor.cells=C;

% Keep only active cells
spin_system.mesh.vor.cells=spin_system.mesh.vor.cells(spin_system.mesh.idx.active);
spin_system.mesh.vor.ncells=numel(spin_system.mesh.vor.cells);

% Voronoi cell area calculation
vor_cell_areas=zeros(spin_system.mesh.vor.ncells,1);
for n=1:spin_system.mesh.vor.ncells
    vor_xcoords=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1);
    vor_ycoords=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2);
    vor_cell_areas(n)=polyarea(vor_xcoords,vor_ycoords);  
end 

% Add weights to mesh structure 
spin_system.mesh.vor.weights=vor_cell_areas;

% Find the maximum number of vertices making up the cell
spin_system.mesh.vor.max_cell_size=max(cellfun(@numel,spin_system.mesh.vor.cells));

end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'idx')
    error('vertex index is missing from spin_system.mesh structure.');
end
end

% A man's rights rest in three boxes: the ballot box,
% the jury box, and the cartridge box.
%
% Frederick Douglass

