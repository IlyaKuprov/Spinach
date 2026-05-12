% Longitudinal spin-phonon coupling producing spin coherence
% modulation and spin-conditioned displacement of a quantised
% vibrational mode. This is a minimal Weyl-algebra version of
% strain-modulated spin Hamiltonians used for NV-centre and
% molecular spin-phonon dynamics.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_dephasing()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V7'};

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
g=2*pi*4e6;

% Longitudinal spin-phonon Hamiltonian
H=omega_v*N+g*Sz*Q;

% Clean up numerical asymmetry
H=(H+H')/2;

% Initial states for coherence and displacement
rho_coh=state(spin_system,{'Lx','BL1'},{1,2});
rho_disp=state(spin_system,{'ZL2','BL1'},{1,2});

% Detection operators
Sx=state(spin_system,'Lx',1);
Qd=state(spin_system,'C',2)+state(spin_system,'A',2);

% Propagate spin coherence and spin-conditioned phonon displacement
traj_s=evolution(spin_system,H,Sx,rho_coh,2e-9,500,'observable');
traj_q=evolution(spin_system,H,Qd,rho_disp,2e-9,500,'observable');

% Validate visible oscillator displacement
if max(abs(real(traj_q)))<0.2
    error('spin-conditioned displacement is not visible.');
end

% Validate visible spin coherence modulation
if (max(real(traj_s))-min(real(traj_s)))<0.02
    error('spin coherence modulation is not visible.');
end

% Plot the coupled dynamics
time_axis=linspace(0,1.0,501);
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,real(traj_s),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$S_X$'); ktitle('spin coherence');
subplot(1,2,2); plot(time_axis,real(traj_q),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, $\mu$s');
kylabel('$a+a^+$'); ktitle('conditional displacement');

end
