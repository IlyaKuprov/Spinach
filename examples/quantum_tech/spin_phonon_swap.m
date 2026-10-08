% Resonant excitation swap between an electron spin and a
% quantised phonon mode. The model is the spin-phonon Jaynes-
% Cummings limit used in mechanical spin-qubit proposals such
% as Rabl et al., Nature Physics 6, 602 (2010).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_swap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V3'};

% Resonant phonon mode in the rotating frame
inter.modes.frqs={[] 0};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=4e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2});
parameters.sweep=1e9;
parameters.npoints=801;

% Trajectory, 'cavity' is the set that keeps spin-mode exchange
traj=device(spin_system,@traject,parameters,'cavity');

% Project out the spin and phonon excitation populations
coil_s=state(spin_system,{'ZL2','E'},{1,2});
coil_v=state(spin_system,{'ZL1','BL2'},{1,2});
pop_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj);
pop_v=cellfun(@(rho)full(hdot(coil_v,rho)),traj);

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
kfigure(); plot(time_axis,real([pop_s; pop_v]),'LineWidth',1.5);
axis tight; ylim([-0.05 1.05]); kgrid; kxlabel('time, ns');
kylabel('excitation population');
ktitle('spin-phonon excitation swap');
klegend({'spin','phonon'},'Location','Best');

end

