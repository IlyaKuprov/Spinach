% Vacuum Rabi swap between a transmon and a microwave cavity
% mode, both represented by truncated bosonic Weyl algebras.
% This is the circuit-QED Jaynes-Cummings limit of Blais et
% al., Rev. Mod. Phys. 93, 025005 (2021).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_cavity_swap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','C3'};

% Resonant transmon-cavity pair in the rotating frame
inter.modes.frqs={0 0};
inter.modes.anharms={-250e6 []};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=20e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,{'BL2','BL1'},{1,2});
parameters.sweep=2e9;
parameters.npoints=301;

% Transmon excitation population through the device context
parameters.coil=state(spin_system,{'BL2','E'},{1,2});
pop_t=device(spin_system,@acquire,parameters,'cavity');

% Cavity excitation population through the device context
parameters.coil=state(spin_system,{'E','BL2'},{1,2});
pop_c=device(spin_system,@acquire,parameters,'cavity');

% Plot the swap dynamics
time_axis=linspace(0,150,301);
kfigure(); plot(time_axis,real([pop_t pop_c]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('transmon-cavity vacuum Rabi swap');
klegend({'transmon','cavity'},'Location','Best');

end

