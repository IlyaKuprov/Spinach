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

% Initial spin excitation and observables
rho=state(spin_system,{'ZL2','BL1'},{1,2});
Ps=state(spin_system,{'ZL2','E'},{1,2});
Pc=state(spin_system,{'ZL1','BL2'},{1,2});

% Propagate the excitation swap
pop_s=evolution(spin_system,H,Ps,rho,1e-9,500,'observable');
pop_c=evolution(spin_system,H,Pc,rho,1e-9,500,'observable');

% Validate visible excitation exchange
if (max(real(pop_c))<0.95)||(min(real(pop_s))>0.05)
    error('vacuum Rabi exchange is not visible.');
end

% Validate population conservation in the active doublet
if max(abs(real(pop_s+pop_c)-1))>1e-6
    error('active-doublet population is not conserved.');
end

% Plot the vacuum Rabi dynamics
time_axis=linspace(0,500,501);
kfigure(); plot(time_axis,real([pop_s pop_c]),'LineWidth',1.5);
axis tight; ylim([-0.05 1.05]); kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('spin-cavity vacuum Rabi oscillation');
klegend({'spin','cavity'},'Location','Best');

end
