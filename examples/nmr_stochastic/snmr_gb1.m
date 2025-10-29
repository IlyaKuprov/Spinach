% A Primas-style stochastic NMR experiment on GB1 protein. The 
% calculation requires a terabyte of RAM and NVidia A100 GPU.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function snmr_gb1()

% Protein data import
options.pdb_mol=1;
options.noshift='delete';
options.select='all';
[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options);

% Magnet field
sys.magnet=18.79;

% Tolerances
sys.tols.inter_cutoff=2.0;
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4; bas.space_level=3;

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='kite';
inter.equilibrium='IME';
inter.tau_c={5e-9};
inter.temperature=298;

% Use GPU arithmetic
% sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the Hamiltonian
H0=hamiltonian(assume(spin_system,'nmr'));

% Get the relaxation superoperator
R=relaxation(spin_system);

% Save the workspace
save('gb1_workspace.mat','-v7.3','-nocompression');

% Control operators
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H'); 
Cx=operator(spin_system,'Lx','13C'); 
Cy=operator(spin_system,'Ly','13C'); 
Nx=operator(spin_system,'Lx','15N');
Ny=operator(spin_system,'Ly','15N');

% Observables
coils=[state(spin_system,'Lx','1H')  state(spin_system,'Ly','1H')  ...
       state(spin_system,'Lx','13C') state(spin_system,'Ly','13C') ...
       state(spin_system,'Lx','15N') state(spin_system,'Ly','15N')];

% Initial condition
H_lab=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H_lab);

% Re-save the workspace
save('gb1_workspace.mat','-v7.3','-nocompression');

% Stochastic process parameters
dt=1e-5;            % seconds 
sig_omega=2*pi*100; % sigma=100 Hz
nsteps=1e6;

% Control noise track
cHx=sig_omega*randn(nsteps,1);
cHy=sig_omega*randn(nsteps,1);
cCx=sig_omega*randn(nsteps,1);
cCy=sig_omega*randn(nsteps,1);
cNx=sig_omega*randn(nsteps,1);
cNy=sig_omega*randn(nsteps,1);

% GPU uploads
L0=gpuArray(H0+1i*R); 
Hx=gpuArray(Hx); Hy=gpuArray(Hy);
Cx=gpuArray(Cx); Cy=gpuArray(Cy);
Nx=gpuArray(Nx); Ny=gpuArray(Ny);
rho=gpuArray(rho); coils=gpuArray(coils);

% Trajectory calculation
report(spin_system,'computing trajectory...'); 
fids=zeros(6,nsteps); tic;
for n=1:nsteps

    % Get observables
    fids(:,n)=gather(coils'*rho);

    % Make evolution generator
    G=L0+cHx(n)*Hx+cHy(n)*Hy+...
         cCx(n)*Cx+cCy(n)*Cy+...
         cNx(n)*Nx+cNy(n)*Ny;

    % Take a time step
    rho=step(spin_system,G,rho,dt);

    % Save checkpoints
    if mod(n,1000)==0
        save('gb1_workspace.mat','-v7.3','-nocompression');
        disp([num2str(n) ' trajectory points done.']);
    end

end

% Performance report
disp(['steps per second: ' num2str(nsteps/toc)]);

% Re-save the workspace
save('gb1_workspace.mat','-v7.3','-nocompression');

% Plotting - control sequences
figure(); scale_figure([2.0 3.0])
time_axis=linspace(0,nsteps*dt,nsteps);
subplot(6,2,1); plot(time_axis,cHx/(2*pi)); ktitle('controls');
axis tight; kgrid; klegend({'Hx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');
subplot(6,2,3); plot(time_axis,cHy/(2*pi)); 
axis tight; kgrid; klegend({'Hy'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');
subplot(6,2,5); plot(time_axis,cCx/(2*pi)); 
axis tight; kgrid; klegend({'Cx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');
subplot(6,2,7); plot(time_axis,cCy/(2*pi)); 
axis tight; kgrid; klegend({'Cy'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');
subplot(6,2,9); plot(time_axis,cNx/(2*pi)); 
axis tight; kgrid; klegend({'Nx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');
subplot(6,2,11); plot(time_axis,cNy/(2*pi)); 
axis tight; kgrid; klegend({'Ny'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('nut. freq., Hz');

% Plotting - trajectories
time_axis=linspace(0,nsteps*dt,nsteps);
subplot(6,2,2); plot(time_axis,fids(1,:)); ktitle('observables');
axis tight; kgrid; klegend({'Hx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');
subplot(6,2,4); plot(time_axis,fids(2,:));
axis tight; kgrid; klegend({'Hy'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');
subplot(6,2,6); plot(time_axis,fids(3,:)); 
axis tight; kgrid; klegend({'Cx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');
subplot(6,2,8); plot(time_axis,fids(4,:));
axis tight; kgrid; klegend({'Cy'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');
subplot(6,2,10); plot(time_axis,fids(5,:)); 
axis tight; kgrid; klegend({'Nx'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');
subplot(6,2,12); plot(time_axis,fids(6,:));
axis tight; kgrid; klegend({'Ny'},'Location','NorthEast');
kxlabel('time, seconds'); kylabel('expect. value');

end

