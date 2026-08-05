% Rabi dynamics of a driven four-level transmon in the Duffing
% approximation, including leakage into the second and third
% excited states. Inspired by standard circuit-QED transmon
% treatments, e.g. Krantz et al., Appl. Phys. Rev. 6, 021318
% (2019).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_rabi_leakage()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T4'};

% Resonantly driven transmon in the rotating frame
inter.modes.frqs={0};
inter.modes.anharms={-250e6};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Resonant drive through the homodecoupling channel of acquire.m
parameters.homodec_oper=(operator(spin_system,'C',1)+...
                         operator(spin_system,'A',1))/2;
parameters.homodec_pwr=25e6;

% Sequence parameters
parameters.rho0=state(spin_system,'BL1',1);
parameters.sweep=1e9;
parameters.npoints=401;

% Level populations through the device context
parameters.coil=state(spin_system,'BL1',1);
p1=device(spin_system,@acquire,parameters,'cavity');
parameters.coil=state(spin_system,'BL2',1);
p2=device(spin_system,@acquire,parameters,'cavity');
parameters.coil=state(spin_system,'BL3',1);
p3=device(spin_system,@acquire,parameters,'cavity');
parameters.coil=state(spin_system,'BL4',1);
p4=device(spin_system,@acquire,parameters,'cavity');

% Plot the leakage dynamics
time_axis=linspace(0,400,401);
kfigure(); plot(time_axis,real([p1 p2 p3 p4]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('level population');
ktitle('transmon Rabi dynamics with leakage');
klegend({'L1','L2','L3','L4'},'Location','Best');

end

