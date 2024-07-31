% Inversion-recovery NOE effect spectrum on strychnine, with the rightmost
% proton signal inverted and a pulse-acquire experiment performed after a
% 500 ms mixing time.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function noe_strychnine()

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

% Build the relaxation superoperator
R=relaxation(spin_system);

% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Start in a state with one spin inverted
Lz9=state(spin_system,{'Lz'},{9});
rho=rho_eq-2*Lz9*(Lz9'*rho_eq)/norm(Lz9);

% Evolve for 500 ms under the relaxation superoperator
rho=evolution(spin_system,1i*R,[],rho,0.5,1,'final');

% Subtract the non-perturbed state
rho=rho-rho_eq;

% Set up a pulse-acquire sequence with the resulting state set as initial
parameters.spins={'1H'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','1H');
parameters.pulse_op=(operator(spin_system,'L+','1H')-...
                     operator(spin_system,'L-','1H'))/2i;
parameters.pulse_angle=pi/2;
parameters.decouple={};
parameters.offset=2800;
parameters.sweep=6500;
parameters.npoints=8192;
parameters.zerofill=65536;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

