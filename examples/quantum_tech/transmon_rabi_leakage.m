% Rabi dynamics of a driven four-level transmon in the Duffing
% approximation, including leakage into the second and third
% excited states. The resonant drive is a part of the rotating
% frame Hamiltonian, and all four level populations come from
% a single trajectory.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_rabi_leakage()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T4'};

% Resonantly driven transmon in the rotating frame
inter.modes.frqs={0};
inter.modes.anharms={-250e6};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Resonant drive on the transmon quadrature
Cr=operator(spin_system,'C',1); 
An=operator(spin_system,'A',1);
H=H+2*pi*25e6*(Cr+An)/2;

% Initial condition
rho=state(spin_system,'BL1',1);

% Trajectory of the driven transmon
traj=evolution(spin_system,H,[],rho,1e-9,400,'trajectory');

% Level populations
pops=zeros(4,401);
for n=1:4
    coil=state(spin_system,['BL' int2str(n)],1);
    pops(n,:)=real(cellfun(@(rho)full(hdot(coil,rho)),traj));
end

% Plot the leakage dynamics
time_axis=linspace(0,400,401);
kfigure(); plot(time_axis,pops','LineWidth',1.5);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('level population');
ktitle('transmon Rabi dynamics with leakage');
klegend({'L1','L2','L3','L4'},'Location','Best');

end

