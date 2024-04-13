% Inversion-recovery NOE effect spectrum on a simple four-spin system, with
% the rightmost proton signal inverted and a pulse-acquire experiment per-
% formed after a very long (five seconds) mixing time. Sequential NOE hops
% with alternating signs are clearly visible in the result.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function noe_four_spin()

% Set the spin system
sys.isotopes={'1H','1H','1H','1H'};
inter.zeeman.scalar={1.0 2.0 3.0 4.0};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,3}=10;
inter.coupling.scalar{3,4}=10;
inter.coordinates={[0.0 0.0 0.0]; [0.0 0.0 2.0];
                   [0.0 0.0 4.0]; [0.0 0.0 6.0]};
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='kite';
inter.temperature=298;
inter.tau_c={200e-12};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build the relaxation superoperator
R=relaxation(spin_system);

% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Start in a state with one spin inverted
Lz1=state(spin_system,{'Lz'},{1});
rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2;

% Evolve for five seconds under the relaxation superoperator
rho=evolution(spin_system,1i*R,[],rho,5.0,1,'final');

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
parameters.offset=1400;
parameters.sweep=4500;
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
scale_figure([2.0 0.85]); xlim([0.8 4.2]);
kylabel('intensity, a.u.');

end

