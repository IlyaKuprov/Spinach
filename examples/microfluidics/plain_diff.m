% Simple diffusion simulation without spin dynamics. Longitudinal 
% magnetisation is tracked as a function of time.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk

function plain_diff()

% One proton
sys.magnet=14.1;
sys.isotopes={'1H'};

% Chemical shift (water)
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic switches
sys.disable={'trajlevel'};
sys.enable={'gpu'};

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

% Initial condition: Lz in one cell in the middle
parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1);
parameters.rho0_ph{1}(1230)=0.5;
parameters.rho0_st{1}=state(spin_system,'Lz','1H');

% Detection state: Lz in all cells
parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1);
parameters.coil_st{1}=state(spin_system,'Lz','1H');

% Sequence and timing parameters
parameters.spins={'1H'};
parameters.offset=0;
parameters.dt=160; 
parameters.npoints=500;

% Set assumptions
spin_system=assume(spin_system,'nmr');

% Same Hamiltonian everywhere
H=hamiltonian(spin_system);
H=frqoffset(spin_system,H,parameters);
parameters.H_op={H}; 
parameters.H_ph={ones(2659,1)};

% Same relaxation everywhere
parameters.R_op={relaxation(spin_system)}; 
parameters.R_ph={ones(2659,1)};

% Just diffusion
spin_system.mesh.u=0*spin_system.mesh.u;
spin_system.mesh.v=0*spin_system.mesh.v;
parameters.diff=1e-7;

% Drainage in the distal pipe
drainage=zeros(2659,1); 
drainage(2600:end)=-0.01;
parameters.K_op={speye([4 4])}; 
parameters.K_ph={drainage};

% Get the trajectory of simple flow
traj=meshflow(spin_system,@simple_flow,parameters);

% Extract the observable quantity
coil=state(spin_system,'Lz','1H');
traj=fpl2phan(traj(:),coil,[2659 parameters.npoints]);

% Make a figure
figure(); scale_figure([1.5 1.5]); 
camproj('perspective'); view(-20,15); axis vis3d;

% Set Z axis extents
spin_system.mesh.zext=[-0.02 0.02];

% Run through trajectory
for n=1:size(traj,2)

    % Do not steal focus
    set(groot,'CurrentFigure',1); cla;

    % Get the concentrations
    conc=full(real(traj(:,n)));
    
    % Update the figure
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc);
    zlim(spin_system.mesh.zext);
    set(gca,'DataAspectRatio',[1 1 0.05]);
    camorbit(0.5,0); drawnow();
        
end 

end

