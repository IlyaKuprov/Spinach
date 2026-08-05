% An exchange-coupled two-electron system with the electrons
% having independent Jaynes-Cummings couplings to the same
% mode of an electromagnetic cavity. A time-domain simulati-
% on starting with transverse spin magnetisation and empty
% cavity mode. Detected on the Lx operator of the spin and
% magnetic field operator of the cavity mode.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function jaynes_cummings_c()

% Magnet field
sys.magnet=0.33;

% System
sys.isotopes={'E','E','C5'};

% Exchange coupling between the electrons
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{1,2}=5e6;

% Cavity resonant with the electrons
inter.modes.frqs={[] [] -sys.magnet*spin('E')/(2*pi)};
inter.modes.exchange=cell(3,3);
inter.modes.exchange{1,3}=2.828e6;
inter.modes.exchange{2,3}=2.728e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.offset=5e6;
parameters.rho0=state(spin_system,{'Lx','BL1'},{1,3})+...
                state(spin_system,{'Lx','BL1'},{2,3});
parameters.sweep=1e8;
parameters.npoints=251;

% Trajectory through the device context
traj=device(spin_system,@traject,parameters,'cavity');

% Detection state, spins
S=state(spin_system,{'Lx'},{1})+...
  state(spin_system,{'Lx'},{2});

% Detection state, cavity mode
B=(state(spin_system,{'C'},{3})-...
   state(spin_system,{'A'},{3}))/2i;

% Project out the observables
traj_s=S'*traj; traj_c=B'*traj;

% Plot the observables
time_axis=linspace(0,2.5,251); % us
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,real(traj_s));
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$S_{X}$ of the spin');
ktitle('spin dynamics');
subplot(1,2,2); plot(time_axis,real(traj_c));
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$(a^{+}-a)/2i$ of the mode');
ktitle('cavity mode dynamics');

end

