% Simple flow simulation with no dynamics in the spin subspace: 
% longitudinal magnetisation is tracked as a function of time af-
% ter injection into the flow field imported from COMSOL with a
% diffusion term also present. The tail of the pipe has drainage
% terms set up using a kinetics superoperator phantom.
%
% a.acharya@soton.ac.uk
% sylwia.ostrowska@kit.edu
% marcel.utz@kit.edu
% ilya.kuprov@weizmann.ac.il

function plain_flow()

% Import hydrodynamics information
comsol.mesh_file='chip_mesh.txt';
comsol.velo_file='chip_velo.txt';
comsol.crop={[286.8 287.5],[576.0 579.0]};
comsol.inactivate=[9 10 19 30 20 25 14 13   ...
                   3372 3373 3380 3381 3382 ...
                   3386 3169 3185 3201 3054 ...
                   3077 3055 3053 3078 3186 ...
                   3168 875 899 897 877 876 ...
                   860 858 885 859 883];
mesh=comsol_import(comsol);

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

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system.mesh=mesh;

% Initial condition: Lz in a few cells
parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1);
parameters.rho0_ph{1}(140:160)=2;
parameters.rho0_st{1}=state(spin_system,'Lz','1H');

% Detection state: Lz in all cells
parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1);
parameters.coil_st{1}=state(spin_system,'Lz','1H');

% Sequence and timing parameters
parameters.spins={'1H'};
parameters.offset=0;
parameters.dt=50; 
parameters.npoints=200;

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

% Diffusion coefficient
parameters.diff=1e-7; % m^2/s

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
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.05]);
    camorbit(0.5,0); drawnow();
        
end 

end

