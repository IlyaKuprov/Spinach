% Nuclear overhauser effect in a heteronuclear two-spin system in
% the short correlation time case.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function noe_two_spin_het()

% Set the spin system
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={0.0 0.0};
inter.coordinates={[0.00 0.00 0.00];
                   [0.00 0.00 1.03]};

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='kite';
inter.temperature=298;
inter.tau_c={100e-12};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build the relaxation superoperator
R=relaxation(spin_system);

% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Start in a state with proton spin inverted
Lz1=state(spin_system,{'Lz'},{1});
rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2;

% Compute the evolution trajectory
coil=[state(spin_system,{'Lz'},{1})...
      state(spin_system,{'Lz'},{2})];
answer=evolution(spin_system,1i*R,coil,rho,1e-2,400,'multichannel');

% Plot the longitudinal magnetization
figure(); plot(linspace(0,4,401),answer);
kxlabel('time, seconds'); kgrid;
kylabel('$S_{\rm{Z}}$ expectation value');
legend({'Proton','Carbon'},'Location','southeast');

end

