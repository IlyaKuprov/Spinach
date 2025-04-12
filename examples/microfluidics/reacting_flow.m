% Flow in the absence of spin dynamics, but presence of two 
% unidirectional second-order chemical reactions.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function reacting_flow()

% Import Diels-Alder cycloaddition
[sys,inter,bas]=dac_reaction();

% Magnet field
sys.magnet=14.1;

% This needs a GPU
sys.enable={'greedy','gpu'};

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

% Rate constants, mol/(L*s)
k1=2.0;  % towards exo  
k2=1.0;  % towards endo

% Cycloaddition reaction generator, including solvent
K=@(x)([-k1*x(2)-k2*x(2)  0                0      0     0;      
         0               -k1*x(1)-k2*x(1)  0      0     0;           
         0                k1*x(1)          0      0     0;
         0                k2*x(1)          0      0     0; 
         0                0                0      0     0]);  

% Strong diffusion
parameters.diff=1e-7;

% Timing parameters
dt=20; npoints=280;

% Get diffusion and flow generator
GF=flow_gen(spin_system,parameters);

% Trajectory preallocation and the initial state
traj=zeros(5,spin_system.mesh.vor.ncells,npoints+1);
traj(1,1240,1)=0.50; traj(2,1246,1)=0.25;

% Time evolution loop
for n=1:npoints

    % Keep the user informed
    report(spin_system,['chemistry time step ' int2str(n) ...
                        '/' int2str(npoints)]);

    % Build a kinetics generator in each cell
    GK=zeros([5 5 spin_system.mesh.vor.ncells],'like',1i);
    parfor k=1:spin_system.mesh.vor.ncells
        GK(:,:,k)=K(traj(:,k,n));
    end

    % Assemble the evolution generator
    G=1i*spblkdiag(GK)+1i*kron(GF,speye(5));

    % Take the time step
    x_curr=traj(:,:,n); x_curr=x_curr(:);
    x_next=step(spin_system,G,x_curr,dt);
    traj(:,:,n+1)=reshape(x_next,[5 spin_system.mesh.vor.ncells]);

end

% Make a figure
figure(); scale_figure([2.5 2.5]);
ksubplot(2,2,1); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,2); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,3); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,4); camproj('perspective'); view(-20,15); axis vis3d;

% Set Z axis extents
spin_system.mesh.zext=[-0.005 0.01];

% Run through trajectory
for n=1:size(traj,3)

    % First reactant
    ksubplot(2,2,1); 
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(1,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0); ktitle('cyclopentadiene');

    % Second reactant
    ksubplot(2,2,2); 
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(2,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0); ktitle('acrylonitrile');

    % First product
    ksubplot(2,2,3); 
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(3,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0); ktitle('exo-NBCN');

    % Second product
    ksubplot(2,2,4); 
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(4,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0); ktitle('endo-NBCN');

    % Draw the frame
    drawnow();
        
end
   
end

