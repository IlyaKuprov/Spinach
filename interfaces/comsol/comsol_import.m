% COMSOL 2D mesh data import, cropping and preprocessing for
% Spinach. Syntax:
%
%                mesh=comsol_import(comsol)
%
% Parameters:
%
%    comsol.mesh_file - name of an ASCII file with 
%                       vertex coordinates and edge 
%                       index produced by COMSOL
%
%    comsol.velo_file - name of an ASCII file with
%                       vertex-centred flow veloci-
%                       ties produced by COMSOL
%
%    comsol.crop      - {[xmin xmax],[ymin ymax]} 
%                       region of the mesh to retain
%
%    comsol.inactivate - a vector with mesh vertex
%                        indices to deactivate
%
% Outputs:
%
%    mesh             – Spinach mesh object with
%
%                      ▸ geometry (.x, .y, .idx)
%                      ▸ flow velocities (.u, .v)
%                      ▸ Voronoi tessellation (.vor)
%                      ▸ fast-plot auxiliaries (.plot)
%
% Notes: internally this routine is just a convenience wrapper
%        that calls the following functions
%
%          ▸ comsol_mesh()      – mesh import
%          ▸ comsol_velo()      – velocity import
%          ▸ mesh_crop()        – region trimming
%          ▸ mesh_vorn()        – Voronoi tessellation
%          ▸ mesh_preplot()     – plotting accelerators
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=comsol_import.m>

function mesh=comsol_import(comsol)

% Check consistency
grumble(comsol);

% Import the mesh
mesh=comsol_mesh(comsol.mesh_file);

% Import the velocities
mesh=comsol_velo(mesh,comsol.velo_file);

% Crop to the region of interest
mesh=mesh_crop(mesh,comsol.crop); 

% Inactivate user-specified vertices
mesh=mesh_inact(mesh,comsol.inactivate);   

% Run Voronoi tessellation
mesh=mesh_vorn(mesh);

% Run graphical output preprocessing
mesh=mesh_preplot(mesh);                           

end

% Consistency enforcement
function grumble(comsol)
if ~isstruct(comsol), error('comsol must be a structure.'); end
if (~isfield(comsol,'mesh_file'))||(~ischar(comsol.mesh_file))
    error('comsol.mesh_file must be a character string.');
end
if (~isfield(comsol,'velo_file'))||(~ischar(comsol.velo_file))
    error('comsol.velo_file must be a character string.');
end
if ~isfield(comsol,'crop'), error('comsol.crop must be specified.'); end
if (~iscell(comsol.crop))||(numel(comsol.crop)~=2)
    error('comsol.crop must be {[xmin xmax],[ymin ymax]}.');
end
if (~isnumeric(comsol.crop{1}))||(~isreal(comsol.crop{1}))||                ...
   (numel(comsol.crop{1}) ~= 2)||(comsol.crop{1}(1) >= comsol.crop{1}(2))|| ...
   (~isnumeric(comsol.crop{2}))||(~isreal(comsol.crop{2}))||                ...
   (numel(comsol.crop{2}) ~= 2)||(comsol.crop{2}(1) >= comsol.crop{2}(2))
    error('comsol.crop must be {[xmin xmax],[ymin ymax]}.');
end
if ~isfield(comsol,'inactivate'), error('comsol.inactivate must be specified.'); end
if (~isnumeric(comsol.inactivate))||(~isvector(comsol.inactivate))|| ...
    any(comsol.inactivate(:)<1)||any(mod(comsol.inactivate(:),1)~=0)
    error('comsol.inactivate must be a vector of positive integers.');
end
end

% "Publish, perish... no, Sir! Publish 
%  rubbish and flourish."
%
% Overheard at a quantum 
% technology conference

