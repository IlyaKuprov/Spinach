% Vacuum Rabi oscillation between an electron spin and a micro-
% wave cavity mode in the Jaynes-Cummings approximation. This
% is the one-spin limit of the spin-ensemble cavity experiments
% of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105,
% 140501 and 140502 (2010).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_cavity_vacuum_rabi()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','C3'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Coupling parameters
delta=0;
g=2*pi*8e6;

% Jaynes-Cummings Hamiltonian
H=delta*operator(spin_system,'Lz',1)+...
  g*(operator(spin_system,{'L+','A'},{1,2})+...
     operator(spin_system,{'L-','C'},{1,2}));

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial state and observables
rho=state(spin_system,{'ZL1','BL1'},{1,2});
Ps=state(spin_system,{'ZL1','E'},{1,2});
Pc=state(spin_system,{'E','BL2'},{1,2});

% Propagate the excitation swap
pop_s=evolution(spin_system,H,Ps,rho,1e-9,500,'observable');
pop_c=evolution(spin_system,H,Pc,rho,1e-9,500,'observable');

% Plot the vacuum Rabi dynamics
time_axis=linspace(0,500,501);
kfigure(); plot(time_axis,real([pop_s pop_c]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('spin-cavity vacuum Rabi oscillation');
klegend({'spin','cavity'},'Location','East');

end

