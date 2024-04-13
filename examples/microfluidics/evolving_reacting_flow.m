% Flow in the presence of spin dynamics and unidirectional 
% first-order chemical reaction. The tail of the pipe has
% drainage terms set up using a relaxation superoperator
% phantom.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk

function evolving_reacting_flow()

% One proton either side
sys.magnet=14.1;
sys.isotopes={'1H','1H'};

% Chemical shifts
inter.zeeman.scalar={0.0235 -0.0110};

% Irreversible reaction
inter.chem.parts={1,2};
inter.chem.rates=[-1e-3  0.0 ; 
                   1e-3  0.0 ];
inter.chem.concs=[1.0  0.0];

% T1/T2 relaxation
inter.relaxation={'t1_t2'};
inter.r1_rates={0.01 0.01};
inter.r2_rates={0.01 0.01};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1;

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

% Initial condition: some Lz at the entrance
Y=spin_system.mesh.y(spin_system.mesh.idx.active);
parameters.rho0_st{1}=state(spin_system,{'Lz'},{1},'chem');
parameters.rho0_ph{1}=0.1*double(Y>578.4);

% No detection stage inside the sequence
parameters.coil_ph={}; parameters.coil_st={};

% Sequence and timing parameters
parameters.spins={'1H'};
parameters.offset=0;
parameters.dt=1; 
parameters.npoints=500;

% Set assumptions
spin_system=assume(spin_system,'nmr');

% Same Hamiltonian everywhere
H=hamiltonian(spin_system);
H=frqoffset(spin_system,H,parameters);
H=H+1e-3*operator(spin_system,'Lx','1H');
parameters.H_op={H}; 
parameters.H_ph={ones(2659,1)};

% Relaxation in the distal pipe
parameters.R_op={relaxation(spin_system)}; 
parameters.R_ph={double(Y<576.3)};

% Kinetics active in a stripe
parameters.K_op={kinetics(spin_system)}; 
parameters.K_ph={double((Y>577.5)&(Y<578))};

% Get the trajectory of simple flow
traj_rho=meshflow(spin_system,@simple_flow,parameters);

% Get trajectory of concentration only, as dictated by kinetics 
traj_conc=meshflow(spin_system,@just_flow,parameters);

% Detection states
P1x=state(spin_system,{'Lx'},{1});
P1y=state(spin_system,{'Ly'},{1});
P1z=state(spin_system,{'Lz'},{1});
P2x=state(spin_system,{'Lx'},{2});
P2y=state(spin_system,{'Ly'},{2});
P2z=state(spin_system,{'Lz'},{2});

% Extract concentration trajectory 
C1=fpl2phan(traj_conc(:),P1z,[2659 parameters.npoints]);
C2=fpl2phan(traj_conc(:),P2z,[2659 parameters.npoints]);

% Extract the observable quantities
P1x=fpl2phan(traj_rho(:),P1x,[2659 parameters.npoints]);
P1y=fpl2phan(traj_rho(:),P1y,[2659 parameters.npoints]);
P1z=fpl2phan(traj_rho(:),P1z,[2659 parameters.npoints]);
P2x=fpl2phan(traj_rho(:),P2x,[2659 parameters.npoints]);
P2y=fpl2phan(traj_rho(:),P2y,[2659 parameters.npoints]);
P2z=fpl2phan(traj_rho(:),P2z,[2659 parameters.npoints]);

% Extract phase information
phase_a=atan2(real(P1y),real(P1x));
phase_b=atan2(real(P2y),real(P2x));

% Extract amplitudes 
amp_a=full(real(sqrt(P1x.^2+P1y.^2+P1z.^2)));
amp_b=full(real(sqrt(P2x.^2+P2y.^2+P2z.^2)));

% Make a figure
figure(); scale_figure([4.50 2.25]);  
subplot(1,2,1); camproj('perspective'); view(-20,15); axis vis3d;
subplot(1,2,2); camproj('perspective'); view(-20,15); axis vis3d;

% Set floor and ceiling
spin_system.mesh.zext=[-0.02 0.02];

% Run through trajectory
for n=1:size(C1,2)

    % Do not steal focus
    set(groot,'CurrentFigure',1); 
    
    % Plot substance A
    subplot(1,2,1); cla;
    conc=full(real(C1(:,n)));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc,[phase_a(:,n) amp_a(:,n)]);
    zlim(spin_system.mesh.zext);
    set(gca,'DataAspectRatio',[1 1 0.05])
    camorbit(0.5,0); ktitle('Lz, substance A');

    % Plot substance B
    subplot(1,2,2); cla;
    conc=full(real(C2(:,n)));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc,[phase_b(:,n) amp_b(:,n)]);
    zlim(spin_system.mesh.zext);
    set(gca,'DataAspectRatio',[1 1 0.05])
    camorbit(0.5,0); ktitle('Lz, Substance B');
    drawnow();
  
end 

end

