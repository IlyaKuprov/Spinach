% Powder spectrum of Gd(III) with ZFS up to 3rd spherical rank
% using the giant spin Hamiltonian formalism in a sweepable 400 
% MHz NMR magnet and microwaves at 263.2 GHz.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function lanthanide_powder()

% Spin system properties
sys.isotopes={'E8'};
inter.zeeman.scalar={1.9918};
inter.giant.coeff={{[0 0 0],[0 0 -4.65e8 0 0],[1e7 0 0 2e7 0 0 1e7]}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0]}};

% Field sweep
sys.magnet=1;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E8'};
parameters.grid='rep_2ang_100pts_sph';
parameters.mw_freq=263.2e9; 
parameters.fwhm=2e-4;
parameters.int_tol=10.0;
parameters.tm_tol=0.1;
parameters.window=[9.32 9.56];
parameters.npoints=4096;
parameters.rspt_order=Inf;

% Run the simulation in the high-T approximation
parameters.rho0=-state(spin_system,'Lz','E8');
[spec,parameters]=fieldsweep(spin_system,parameters);

% Plotting
kfigure(); plot(parameters.b_axis,spec);
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
xlim tight; ylim padded; kgrid;

end

