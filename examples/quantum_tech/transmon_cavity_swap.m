% Vacuum Rabi swap between a transmon and a microwave cavity
% mode, both represented by truncated bosonic Weyl algebras.
% This is the circuit-QED Jaynes-Cummings limit of Blais et
% al., Phys. Rev. A 69, 062320 (2004).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_cavity_swap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','C3'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Transmon and cavity operators
Nt=operator(spin_system,'N',1);
Kt=operator(spin_system,'CCAA',1);

% Coupling parameters
anharm=-2*pi*250e6;
delta=0;
g=2*pi*20e6;

% Transmon-cavity exchange Hamiltonian
H=delta*Nt+(anharm/2)*Kt+...
  g*(operator(spin_system,{'C','A'},{1,2})+...
     operator(spin_system,{'A','C'},{1,2}));

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial state and observables
rho=state(spin_system,{'BL2','BL1'},{1,2});
Pt=state(spin_system,{'BL2','E'},{1,2});
Pc=state(spin_system,{'E','BL2'},{1,2});

% Propagate the excitation swap
pop_t=evolution(spin_system,H,Pt,rho,0.5e-9,300,'observable');
pop_c=evolution(spin_system,H,Pc,rho,0.5e-9,300,'observable');

% Plot the swap dynamics
time_axis=linspace(0,150,301);
kfigure(); plot(time_axis,real([pop_t pop_c]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('transmon-cavity vacuum Rabi swap');
klegend({'transmon','cavity'},'Location','Best');

end

