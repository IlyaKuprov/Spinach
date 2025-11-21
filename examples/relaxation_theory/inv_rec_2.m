% An example of inversion recovery experiment simulation 
% for a strychnine spin system.
%
% Calculation time: minutes.
%
% mariagrazia.concilio@sjtu.edu.cn

function inv_rec_2()

% Read the spin system properties 
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=14.1;

% Disable Krylov propagation
sys.disable={'krylov'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='kite';
inter.temperature=298;
inter.tau_c={200e-12};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Aquisition parameters
parameters.sweep=6500; 
parameters.npoints=8192;  
parameters.zerofill=65536;
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.offset=2800;   % Hz

% Set up different recovery delays
mixing_time=[0.01; 0.1; 0.5; 1; 5; 10]; % s

% Lab frame Hamiltonian
H=hamiltonian(assume(spin_system,'labframe'),'left');

% Initial state - thermal equilibrium
parameters.rho_eq=equilibrium(spin_system,H);

% Detection state
parameters.coil=state(spin_system,'L+','1H');

% Rotating frame Liouvillian
L=hamiltonian(assume(spin_system,'nmr'))+...
  1i*relaxation(spin_system);

% Pulse operator
Ly=operator(spin_system,'Ly','1H');

% Get figure going
kfigure(); scale_figure([2.5 1.8]);

% Set up loop
for n=1:numel(mixing_time)
    
    % Run different recovery delays
    parameters.tau=mixing_time(n);

    % Apply the 180 degree inversion pulse
    rho=step(spin_system,Ly,parameters.rho_eq,pi);

    % Evolution 
    rho=evolution(spin_system,L,[],rho,parameters.tau,1,'final');

    % Apply the 90 degree pulse
    rho=step(spin_system,Ly,rho,pi/2);

    % Acquisition
    fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep,parameters.npoints-1,'observable');

    % Apodisation
    fid=apodisation(spin_system,fid,{{'exp',5}});

    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));

    % Plotting
    subplot(2,3,n); plot_1d(spin_system,real(spectrum),parameters);  
    ktitle(['recovery delay = ' num2str(mixing_time(n)), 's']); drawnow();
       
end
 
end

