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

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Zeeman interaction operator, spin
Ez=operator(spin_system,{'Lz'},{1});

% Jaynes-Cummings term
g=2*pi*2.828e6;
H_JC=g*(operator(spin_system,{'L+','A'},{1,2})+...
        operator(spin_system,{'L-','C'},{1,2}));

% Detuning of 5 MHz
delta=2*pi*5e6;

% Rotating frame Hamiltonian
H=delta*Ez+H_JC;

% Initial state - uncorrelated Sx on spin
% and the cavity mode initially empty
rho=state(spin_system,{'Lx','E'  },{1,2})+...
    state(spin_system,{'E' ,'BL1'},{1,2})/2;

% Detection state, spin
S=state(spin_system,{'Lx'},{1});

% Detection state, cavity mode
B=(state(spin_system,{'C'},{2})-...
   state(spin_system,{'A'},{2}))/2i;

% Run observables evolution
traj_s=evolution(spin_system,H,S,rho,10e-9,250,'observable');
traj_c=evolution(spin_system,H,B,rho,10e-9,250,'observable');

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
L1=evolution(spin_system,H,L1,rho,10e-9,250,'observable');
L2=state(spin_system,{'E','BL2'},{1,2});
L2=evolution(spin_system,H,L2,rho,10e-9,250,'observable');
L3=state(spin_system,{'E','BL3'},{1,2});
L3=evolution(spin_system,H,L3,rho,10e-9,250,'observable');
L4=state(spin_system,{'E','BL4'},{1,2});
L4=evolution(spin_system,H,L4,rho,10e-9,250,'observable');
L5=state(spin_system,{'E','BL5'},{1,2});
L5=evolution(spin_system,H,L5,rho,10e-9,250,'observable');
kfigure(); plot(time_axis,real([L1 L2 L3 L4 L5])); 
xlim tight; ylim padded; kxlabel('time, $\mu$s');
kylabel('cavity level population'); kgrid;
klegend({'L1','L2','L3','L4','L5'},'Location','East');

end

