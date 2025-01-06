% Powder averaged W-band field-swept ESR spectrum of Gd(III)
% DOTA complex. Exact diagonalisation is used.
%
% Calculation time: seconds.
%
% corzilius@solidstatednp.com
% ilya.kuprov@weizmann.ac.il

function fieldsweep_gd_dota()

% Isotopes
sys.isotopes={'E8'};

% Magnet field (must be 1)
sys.magnet=1;

% Properties
inter.zeeman.scalar={1.9918};
inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3;
inter.coupling.euler{1,1}=[0 0 0];

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E8'};
parameters.grid=4;
parameters.mw_freq=90e9;
parameters.fwhm=2e-4;
parameters.int_tol=0.01;
parameters.tm_tol=0.1;
parameters.window=[3.05 3.4];
parameters.npoints=4096;
parameters.rspt_order=Inf;

% Run the simulation in the high-T approximation
parameters.rho0=-state(spin_system,'Lz','E8');
[b_axis,spec]=fieldsweep(spin_system,parameters);

% Plotting
figure(); plot(b_axis',spec');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
axis tight; kgrid;

end

