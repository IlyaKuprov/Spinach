% Vacuum Rabi swap between a transmon and a microwave cavity
% mode, both represented by truncated bosonic Weyl algebras.
% This is the circuit-QED Jaynes-Cummings limit of Blais et
% al., Rev. Mod. Phys. 93, 025005 (2021).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_cavity_swap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','C3'};

% Resonant transmon-cavity pair in the rotating frame
inter.modes.frqs={0 0};
inter.modes.anharms={-250e6 []};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=20e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,{'BL2','BL1'},{1,2});
parameters.sweep=2e9;
parameters.npoints=301;

% Trajectory through the device context
traj=device(spin_system,@traject,parameters,'cavity');

% Project out the transmon and cavity excitation populations
coil_t=state(spin_system,{'BL2','E'},{1,2});
coil_c=state(spin_system,{'E','BL2'},{1,2});
pop_t=cellfun(@(rho)full(hdot(coil_t,rho)),traj);
pop_c=cellfun(@(rho)full(hdot(coil_c,rho)),traj);

% Validate visible excitation exchange
if (max(real(pop_c))<0.95)||(min(real(pop_t))>0.05)
    error('transmon-cavity exchange is not visible.');
end

% Validate population conservation in the active doublet
if max(abs(real(pop_t+pop_c)-1))>1e-6
    error('active-doublet population is not conserved.');
end

% Plot the swap dynamics
time_axis=linspace(0,150,301);
kfigure(); plot(time_axis,real([pop_t; pop_c]),'LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('transmon-cavity vacuum Rabi swap');
klegend({'transmon','cavity'},'Location','Best');

end

