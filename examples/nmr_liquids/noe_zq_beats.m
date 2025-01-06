% Zero-quantum beats in the Overhauser effect in a strongly
% coupled two-spin system.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% cthiele@thielelab.de

function noe_zq_beats()

% Set the spin system
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0.0 0.01};           % ppm
inter.coupling.scalar={0.0 3.0; 0.0 0.0}; % Hz
inter.coordinates={[0.00 0.00 0.00]; 
                   [0.00 0.00 2.00]};     % Angstrom
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='secular';
inter.temperature=298;
inter.tau_c={1e-9};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build the Liouvillian
L=hamiltonian(assume(spin_system,'nmr'))+1i*relaxation(spin_system);

% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Start in a state with one spin inverted
Lz1=state(spin_system,{'Lz'},{1});
rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2;

% Compute the evolution trajectory
coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2})];
answer=evolution(spin_system,L,coil,rho,1e-2,1000,'multichannel');

% Plot the longitudinal magnetization
figure(); 
plot(linspace(0,10,1001),real(answer));
kxlabel('time, seconds'); kgrid;
kylabel('$S_{\rm{Z}}$ expectation value');
legend({'Proton A','Proton B'},'Location','southeast');

end

