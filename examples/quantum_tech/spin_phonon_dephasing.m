% Longitudinal spin-phonon coupling producing spin coherence
% modulation by a quantised vibrational mode. This is a minimal
% Weyl-algebra version of strain-modulated spin Hamiltonians
% used for NV-centre and molecular spin-phonon dynamics.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_dephasing()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V5'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Spin and phonon operators
Sz=operator(spin_system,'Lz',1);
N=operator(spin_system,'N',2);
Q=operator(spin_system,'C',2)+operator(spin_system,'A',2);

% Longitudinal spin-phonon Hamiltonian parameters
omega_v=2*pi*20e6;
g=2*pi*3e6;

% Longitudinal spin-phonon Hamiltonian
H=omega_v*N+g*Sz*Q;

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial transverse spin state with the phonon in vacuum
rho=state(spin_system,{'Lx','BL1'},{1,2});
Sx=state(spin_system,'Lx',1);
Qd=state(spin_system,'C',2)+state(spin_system,'A',2);

% Propagate spin coherence and phonon displacement
traj_s=evolution(spin_system,H,Sx,rho,2e-9,500,'observable');
traj_q=evolution(spin_system,H,Qd,rho,2e-9,500,'observable');

% Plot the coupled dynamics
time_axis=linspace(0,1.0,501);
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,real(traj_s),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$S_X$'); ktitle('spin coherence');
subplot(1,2,2); plot(time_axis,real(traj_q),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$a+a^+$'); ktitle('phonon displacement');

end

