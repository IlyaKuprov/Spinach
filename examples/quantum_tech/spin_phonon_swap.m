% Resonant excitation swap between an electron spin and a
% quantised phonon mode. The model is the spin-phonon Jaynes-
% Cummings limit used in mechanical spin-qubit proposals such
% as Rabl et al., Phys. Rev. B 79, 041302 (2009).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_swap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V3'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Coupling parameters
delta=0;
g=2*pi*4e6;

% Spin-phonon exchange Hamiltonian
H=delta*operator(spin_system,'Lz',1)+...
  g*(operator(spin_system,{'L+','A'},{1,2})+...
     operator(spin_system,{'L-','C'},{1,2}));

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial spin excitation and observables
rho=state(spin_system,{'ZL2','BL1'},{1,2});
Ps=state(spin_system,{'ZL2','E'},{1,2});
Pv=state(spin_system,{'ZL1','BL2'},{1,2});

% Propagate the excitation swap
pop_s=evolution(spin_system,H,Ps,rho,1e-9,800,'observable');
pop_v=evolution(spin_system,H,Pv,rho,1e-9,800,'observable');

% Validate visible excitation exchange
if (max(real(pop_v))<0.95)||(min(real(pop_s))>0.05)
    error('spin-phonon exchange is not visible.');
end

% Validate population conservation in the active doublet
if max(abs(real(pop_s+pop_v)-1))>1e-6
    error('active-doublet population is not conserved.');
end

% Plot the swap dynamics
time_axis=linspace(0,800,801);
kfigure(); plot(time_axis,real([pop_s pop_v]),'LineWidth',1.5);
axis tight; ylim([-0.05 1.05]); kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('spin-phonon excitation swap');
klegend({'spin','phonon'},'Location','East');

end
