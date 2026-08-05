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

% Resonant cavity in the rotating frame
inter.modes.frqs={[] 0};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=8e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2});
parameters.sweep=1e9;
parameters.npoints=501;

% Trajectory through the device context
traj=device(spin_system,@traject,parameters,'cavity');

% Project out the spin and cavity excitation populations
coil_s=state(spin_system,{'ZL2','E'},{1,2});
coil_c=state(spin_system,{'ZL1','BL2'},{1,2});
pop_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj);
pop_c=cellfun(@(rho)full(hdot(coil_c,rho)),traj);

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
kfigure(); plot(time_axis,real([pop_s; pop_c]),'LineWidth',1.5);
axis tight; ylim([-0.05 1.05]); kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('spin-cavity vacuum Rabi oscillation');
klegend({'spin','cavity'},'Location','Best');

end

