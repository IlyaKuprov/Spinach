% Jaynes-Cummings coupling between a spin and an electromagnetic
% cavity mode with five population numbers included. A time-dom-
% ain simulation starting with transverse spin magnetisation and
% empty cavity mode. Detected on the Lx operator of the spin and
% magnetic field operator of the cavity mode.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function jaynes_cummings_b()

% Magnet field
sys.magnet=0.33;

% System
sys.isotopes={'E','C5'};

% Cavity resonant with the electron
e_frq=-sys.magnet*spin('E')/(2*pi);
inter.modes.frqs={[] e_frq};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=2.828e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.offset=5e6;
parameters.rho0=state(spin_system,{'Lx','E'  },{1,2})+...
                state(spin_system,{'E' ,'BL1'},{1,2})/2;
parameters.sweep=1e8;
parameters.npoints=251;

% Trajectory through the device context
traj=device(spin_system,@traject,parameters,'cavity');

% Detection state, spin
S=state(spin_system,{'Lx'},{1});

% Detection state, cavity mode
B=(state(spin_system,{'C'},{2})-...
   state(spin_system,{'A'},{2}))/2i;

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

% Cavity energy level population evolution
L1=state(spin_system,{'E','BL1'},{1,2});
L2=state(spin_system,{'E','BL2'},{1,2});
L3=state(spin_system,{'E','BL3'},{1,2});
pop_1=L1'*traj; pop_2=L2'*traj; pop_3=L3'*traj;

% Plot the level populations
kfigure(); plot(time_axis,real([pop_1; pop_2; pop_3]));
kxlabel('time, $\mu$s'); kylabel('population'); kgrid;
axis tight; ktitle('cavity level populations');
klegend({'BL1','BL2','BL3'},'Location','Best');

end

