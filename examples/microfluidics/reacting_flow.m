% Flow in the absence of spin dynamics (Lz is detected), but pre-
% sence of a unidirectional first-order chemical reaction. The 
% tail of the pipe has drainage terms set up using a kinetics 
% superoperator phantom.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk

function reacting_flow()

% One proton on either side
sys.magnet=14.1;
sys.isotopes={'1H','1H'};

% Chemical shifts
inter.zeeman.scalar={0.0 0.0};

% Irreversible reaction
inter.chem.parts={1,2};
inter.chem.rates=[-1e-3  0.0; 
                   1e-3  0.0];
inter.chem.concs=[1.0  0.0];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic switches
sys.disable={'trajlevel'};
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% COMSOL mesh import
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

% Initial condition: Lz of the first spin in a few cells
parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1);
parameters.rho0_ph{1}(140:160)=1;
parameters.rho0_st{1}=state(spin_system,{'Lz'},{1},'chem');

% No detection stage inside the sequence
parameters.coil_ph={}; parameters.coil_st={};

% Sequence and timing parameters
parameters.spins={'1H'};
parameters.offset=0;
parameters.dt=160; 
parameters.npoints=50;

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

% Kinetics A: reaction in a stripe 
Y=spin_system.mesh.y(spin_system.mesh.idx.active);
Ph1=double((Y>577.5)&(Y<578)); K1=kinetics(spin_system);

% Kinetics B: drainage in the distal pipe
Ph2=zeros(2659,1); Ph2(2600:end)=-0.01; 
K2=unit_oper(spin_system);

% Assemble kinetics
parameters.K_op={K1,K2}; 
parameters.K_ph={Ph1,Ph2};

% Get the trajectory of simple flow
traj=meshflow(spin_system,@simple_flow,parameters);

% Extract the observable quantity
coil_a=state(spin_system,{'Lz'},{1});
coil_b=state(spin_system,{'Lz'},{2});
traj_a=fpl2phan(traj(:),coil_a,[2659 parameters.npoints]);
traj_b=fpl2phan(traj(:),coil_b,[2659 parameters.npoints]);

% Make a figure
figure(); scale_figure([4.50 2.25]); 
subplot(1,2,1); camproj('perspective'); view(-20,15); axis vis3d;
subplot(1,2,2); camproj('perspective'); view(-20,15); axis vis3d;

% Set the ceiling
spin_system.mesh.zext=[-0.02 0.02];

% Run through trajectory
for n=1:size(traj,2)

    % Do not steal focus
    set(groot,'CurrentFigure',1); 

    subplot(1,2,1); cla;
    conc=full(real(traj_a(:,n)));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc);
    zlim(spin_system.mesh.zext);
    set(gca,'DataAspectRatio',[1 1 0.05])
    camorbit(0.5,0); ktitle('Substance A');

    subplot(1,2,2); cla;
    conc=full(real(traj_b(:,n)));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc);
    zlim(spin_system.mesh.zext);
    set(gca,'DataAspectRatio',[1 1 0.05])
    camorbit(0.5,0); ktitle('Substance B');
    drawnow();
  
end 

end

