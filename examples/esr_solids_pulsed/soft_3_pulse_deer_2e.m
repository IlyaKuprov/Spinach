% Three-pulse DEER simulation for a two-electron system. Soft 
% pulses are simulated using the Fokker-Planck formalism.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function soft_3_pulse_deer_2e()

% Magnet field
sys.magnet=0.3451805;

% Isotopes
sys.isotopes={'E','E'};

% Zeeman interactions
inter.zeeman.eigs=cell(1,2);
inter.zeeman.euler=cell(1,2);
inter.zeeman.eigs{1}=[2.284 2.123 2.075];
inter.zeeman.euler{1}=[135 90 45]*(pi/180);
inter.zeeman.eigs{2}=[2.035 2.013 1.975];
inter.zeeman.euler{2}=[30 60 120]*(pi/180);

% Coordinates (Angstrom)
inter.coordinates=cell(2,1);
inter.coordinates{1}=[0.00  0.00 0.00];
inter.coordinates{2}=[20.00 0.00 0.00];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.grid='rep_2ang_6400pts_sph';
parameters.method='expm';
parameters.verbose=0;

% EPR parameters
parameters.offset=-4e8;
parameters.sweep=3e9;
parameters.npoints=256;
parameters.zerofill=2048;
parameters.axis_units='GHz-labframe';
parameters.derivative=0;
parameters.invert_axis=1;
parameters.assumptions='deer';

% DEER pulse parameters
parameters.pulse_rnk=[2 2 2];
parameters.pulse_dur=[20e-9 50e-9 40e-9];
parameters.pulse_phi=[pi/2 pi/2 pi/2];
parameters.pulse_pwr=2*pi*[8e6 8e6 8e6];
parameters.pulse_frq=[9.720e9 10.255e9 9.720e9];

% DEER echo timing parameters
parameters.p1_p3_gap=1e-6;
parameters.p2_nsteps=100;
parameters.echo_time=100e-9;
parameters.echo_npts=100;

% Simulation and plotting
deer_3p_soft_diag(spin_system,parameters);

end

