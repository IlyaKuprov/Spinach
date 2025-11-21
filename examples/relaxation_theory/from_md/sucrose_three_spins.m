% One of the calculations reported in the JMR paper with Jim Prestegard:
% a three-spin subsystem from the glucose ring of sucrose - an illustra-
% tion of incorrect viscosity of TIP3P water. 
%
% The experimental value of tau_c for sucrose is around 90 ps (and this
% is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees
% with Redfield theory when tau_c is set to 39 ps in the latter.
%
% Here, numerical relaxation superoperator computed from a long MD tra-
% jectory is compared with the analytical one computed using the isotro-
% pic rotational diffusion approximation.
%
%               https://doi.org/10.1016/j.jmr.2020.106891
%
% Calculation time: minutes, with most of the time spent 
%                   computing MD frame Hamiltonians
%
% jpresteg@uga.edu
% ilya.kuprov@weizmann.ac.il

function sucrose_three_spins()

% Three-spin system
sys.magnet=14.1;
sys.isotopes={'1H','1H','1H'};

% Complete basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Analytical relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={39e-12};
inter.temperature=298;

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set laboratory frame assumptions
spin_system=assume(spin_system,'labframe');

% Get coherent Hamiltonian
H0=hamiltonian(spin_system);

% Load MD trajectory, dims: (XYZ, spins, time)
load('suc_three_spin_traj.mat','traj','dt');

% 50k frames is enough here
traj=traj(:,:,1:50000); ref_frame=1;

% Compute Hamiltonian for each trajectory frame
report(spin_system,'computing MD frame Hamiltonians...');
H1=cell(1,size(traj,3));
parfor n=1:size(traj,3)

    % Localise system object
    spinsys_loc=spin_system;
    
    % Hush up Spinach output
    spinsys_loc.sys.output='hush';
    
    % Pull a trajectory frame
    traj_slice=traj(:,:,n);
    
    % Set current coordinates
    spinsys_loc.inter.coordinates={double(traj_slice(:,1)'); 
                                   double(traj_slice(:,2)'); 
                                   double(traj_slice(:,3)')}; 
    spinsys_loc.inter.pbc={}; spinsys_loc=dipolar(spinsys_loc);

    % Get the anisotropic part of the Hamiltonian
    [~,Q]=hamiltonian(spinsys_loc); H1{n}=orientation(Q,[0,0,0]);
    
end

% Get MD relaxation superoperator
tau_est=inter.tau_c{1}; tic;
R_gce=ngce(spin_system,H0,H1,dt,tau_est,0);
ngce_run_time=toc;

% Load reference frame coordinates
spin_system.inter.coordinates={double(traj(:,1,ref_frame)'); 
                               double(traj(:,2,ref_frame)'); 
                               double(traj(:,3,ref_frame)')};
spin_system.inter.pbc={}; spin_system=dipolar(spin_system);

% Get analytical Redfield matrix
R_red=relaxation(spin_system);

% Plotting
kfigure();
plot(diag(R_gce),diag(R_red),'ro');
kxlabel('Redfield matrix elements, mol. dyn.');
kylabel('Redfield matrix elements, rot. diff.');
edges=[min([diag(R_gce); diag(R_red)]),...
       max([diag(R_gce); diag(R_red)])];
xlim(edges); ylim(edges); 
kgrid; box on; axis square;
disp(['NGCE stage run time, seconds: ' ...
      num2str(ngce_run_time)]);

end

