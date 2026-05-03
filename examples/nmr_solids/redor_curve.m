% REDOR dephasing curve for a simple 13C-15N spin pair using
% Fokker-Planck MAS formalism.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function redor_curve()

% System specification
sys.magnet=14.1;
sys.isotopes={'13C','15N'};

% Interactions
inter.zeeman.scalar={0.0 0.0};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.47]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};
sys.enable={'prop_cache'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% REDOR setup
parameters.spins={'13C','15N'};
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=9;
parameters.grid='leb_2ang_rank_11';
parameters.ncycles=0:2:24;
parameters.rho0=state(spin_system,'Lx','13C');
parameters.coil=state(spin_system,'Lx','13C');
parameters.verbose=0;

% Simulation
curve=singlerot(spin_system,@redor,parameters,'nmr');

% Normalised REDOR difference
redor_diff=real(curve(3,:)./curve(1,:));

% Plotting
kfigure(); plot(parameters.ncycles,redor_diff,'o-'); kgrid;
kxlabel('REDOR evolution time, rotor cycles');
kylabel('$\Delta S/S_0$'); xlim tight;

end

