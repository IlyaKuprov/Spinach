% Import, Voronoi tessellation, and plotting of the 
% hydrodynamic mesh and velocity field from COMSOL. 
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function show_mesh()

% Import hydrodynamics information
comsol.mesh_file='mesh-4ulm.txt';
comsol.velo_file='velocity-field-4ulm.txt';
comsol.crop={[286.8 287.5],[576.0 579.0]};
comsol.inactivate=[9 10 19 30 20 25 14 13   ...
                   3372 3373 3380 3381 3382 ...
                   3386 3169 3185 3201 3054 ...
                   3077 3055 3053 3078 3186 ...
                   3168 875 899 897 877 876 ...
                   860 858 885 859 883];
mesh=comsol_import(comsol);

% No spin system here
spin_system=bootstrap();
spin_system.mesh=mesh;

% Draw the mesh
figure(); mesh_plot(spin_system,2,0); 
xlim([286.88  287.42]); ylim([578.07  578.50]);
klegend({'triangles','rectangles',...
         'tessellation','velocities'},'Location','NorthEast');

end

