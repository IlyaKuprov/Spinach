% Nuclear overhauser effect in a homonuclear two-spin system in
% the long correlation time case.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function noe_two_spin_hom()

% Set the spin system
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0.0 0.0};
inter.coordinates={[0.00 0.00 0.00];
                   [0.00 0.00 2.00]};

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
inter.tau_c={1e-9};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build the relaxation superoperator
R=relaxation(spin_system);

% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Start in a state with one spin inverted
Lz1=state(spin_system,{'Lz'},{1});
rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2;

% Compute the evolution trajectory
coil=[state(spin_system,{'Lz'},{1})...
      state(spin_system,{'Lz'},{2})];
answer=evolution(spin_system,1i*R,coil,rho,1e-2,1000,'multichannel');

% Plot the longitudinal magnetization
kfigure(); plot(linspace(0,10,1001),answer);
kxlabel('time, seconds'); kgrid;
kylabel('$S_{\rm{Z}}$ expectation value');
klegend({'Proton A','Proton B'},'Location','southeast');

end

