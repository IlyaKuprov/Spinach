% Flow in the absence of spin dynamics, but presence of two 
% unidirectional second-order chemical reactions.
%
% Simulation time: seconds.
%
% a.acharya@soton.ac.uk
% sylwia.ostrowska@kit.edu
% marcel.utz@kit.edu
% ilya.kuprov@weizmann.ac.il

function reacting_flow()

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

% No spin system here
spin_system=bootstrap();
spin_system.mesh=mesh;

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
kfigure(); scale_figure([2.5 2.5]);
ksubplot(2,2,1); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,2); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,3); camproj('perspective'); view(-20,15); axis vis3d;
ksubplot(2,2,4); camproj('perspective'); view(-20,15); axis vis3d;

% Run through trajectory
for n=1:size(traj,3)

    % First reactant
    ksubplot(2,2,1); 
    spin_system.mesh.zext=[-0.005 0.01];
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(1,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0.02); ktitle('cyclopentadiene');

    % Second reactant
    ksubplot(2,2,2); 
    spin_system.mesh.zext=[-0.005 0.01];
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(2,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025]);
    camorbit(0.5,0.02); ktitle('acrylonitrile');

    % First product
    ksubplot(2,2,3); 
    spin_system.mesh.zext=[-0.005 0.01]/8;
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(3,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025/8]);
    camorbit(0.5,0.02); ktitle('exo-NBCN');

    % Second product
    ksubplot(2,2,4); 
    spin_system.mesh.zext=[-0.005 0.01]/8;
    set(groot,'CurrentFigure',1); cla;
    conc=squeeze(full(real(traj(4,:,n))));
    mesh_plot(spin_system,0,0);
    conc_plot(spin_system,conc');
    zlim(spin_system.mesh.zext);
    kzlabel('concentration, a.u.');
    set(gca,'DataAspectRatio',[1 1 0.025/8]);
    camorbit(0.5,0.02); ktitle('endo-NBCN');

    % Draw the frame
    drawnow();
        
end
   
end

