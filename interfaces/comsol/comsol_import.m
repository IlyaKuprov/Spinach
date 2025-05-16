% The headline COMSOL data import function.
%
% ilya.kuprov@weizmann.ac.il

function mesh=comsol_import(comsol)

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

% "Publish, perish... no, Sir! Publish 
%  rubbish and flourish."
%
% Overheard at a quantum 
% technology conference

