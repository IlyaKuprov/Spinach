% Rabi dynamics of a driven four-level transmon in the Duffing
% approximation, including leakage into the second and third
% excited states. Inspired by standard circuit-QED transmon
% control models.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_rabi_leakage()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T4'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Transmon operators
Cr=operator(spin_system,'C',1);
An=operator(spin_system,'A',1);
K=operator(spin_system,'CCAA',1);

% Rotating-frame drive parameters
anharm=-2*pi*250e6;
omega_1=2*pi*25e6;

% Driven Duffing Hamiltonian
H=(anharm/2)*K+(omega_1/2)*(Cr+An);

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial state and level projectors
rho=state(spin_system,'BL1',1);
L1=state(spin_system,'BL1',1);
L2=state(spin_system,'BL2',1);
L3=state(spin_system,'BL3',1);
L4=state(spin_system,'BL4',1);

% Propagate level populations
p1=evolution(spin_system,H,L1,rho,1e-9,400,'observable');
p2=evolution(spin_system,H,L2,rho,1e-9,400,'observable');
p3=evolution(spin_system,H,L3,rho,1e-9,400,'observable');
p4=evolution(spin_system,H,L4,rho,1e-9,400,'observable');

% Plot the leakage dynamics
time_axis=linspace(0,400,401);
kfigure(); plot(time_axis,real([p1 p2 p3 p4]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('level population');
ktitle('transmon Rabi dynamics with leakage');
klegend({'L1','L2','L3','L4'},'Location','Best');

end

