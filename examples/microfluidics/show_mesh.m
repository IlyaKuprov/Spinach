% Import, Voronoi tessellation, and plotting of the 
% hydrodynamic mesh and velocity field from COMSOL. 
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk

function show_mesh()

% One proton
sys.magnet=14.1;
sys.isotopes={'1H'};

% Chemical shift (water)
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% COMSOL mesh and velocity field import
spin_system=comsol_mesh(spin_system,'mesh-4ulm.txt');            % Read the mesh
spin_system=comsol_velo(spin_system,'velocity-field-4ulm.txt');  % Read velocities
spin_system=mesh_crop(spin_system,[286.8 287.5],[576.0 579.0]);  % Crop the mesh
spin_system=mesh_inact(spin_system,[9 10 19 30 20 25 14 13   ...
                                    3372 3373 3380 3381 3382 ...
                                    3386 3169 3185 3201 3054 ... % Prune out edge vertices
                                    3077 3055 3053 3078 3186 ...
                                    3168 875 899 897 877 876 ...
                                    860 858 885 859 883]);      
spin_system=mesh_vorn(spin_system);                              % Run Voronoi tessellation
spin_system=mesh_preplot(spin_system);                           % Run output preprocessing

% Draw the mesh
figure(); mesh_plot(spin_system,2,0); 
xlim([286.88  287.42]); ylim([578.07  578.50]);
klegend({'triangles','rectangles',...
         'tessellation','velocities'},'Location','NorthEast');

end

