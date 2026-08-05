% Longitudinal spin-phonon coupling producing spin coherence
% modulation and spin-conditioned displacement of a quantised
% vibrational mode. This is a minimal Weyl-algebra version of
% strain-modulated spin Hamiltonians used for NV-centre and
% molecular spin-phonon dynamics. Both observables come from a
% single trajectory: the drift commutes with the spin projec-
% tion operator, so the coherence and the population sectors
% of the initial condition evolve and are detected without
% mixing.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_dephasing()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V7'};

% Phonon mode with a longitudinal coupling
inter.modes.frqs={[] 20e6};
inter.modes.longitudinal=cell(2,2);
inter.modes.longitudinal{1,2}=4e6*sqrt(2);

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.sweep=5e8;
parameters.npoints=501;
parameters.rho0=state(spin_system,{'Lx','BL1'},{1,2})+...
                state(spin_system,{'ZL2','BL1'},{1,2});

% Trajectory through the device context
traj=device(spin_system,@traject,parameters,'spin-phonon');

% Project out the spin coherence and the conditional displacement
coil_s=state(spin_system,'Lx',1);
coil_q=state(spin_system,'C',2)+state(spin_system,'A',2);
traj_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj);
traj_q=cellfun(@(rho)full(hdot(coil_q,rho)),traj);

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

