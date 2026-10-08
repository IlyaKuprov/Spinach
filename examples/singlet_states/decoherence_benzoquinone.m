% Long-lived spin states in the para-benzoquinone molecule
% (4 protons, 256-dimensional Liouville space). The relaxation
% superoperator accounts for every dipolar coupling and every
% CSA tensor in the system.
%
% Calculation time: seconds
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function decoherence_benzoquinone()

% Read the spin system (coordinates, chemical shifts,
% J-couplings and CSAs) from a vacuum DFT calculation
[sys,inter]=g2spinach(gparse('../standard_systems/benzoquinone.log'),...
                                               {{'H','1H'}},31.8,[]);
% Set magnet field to 1.0 Tesla
sys.magnet=1.0;

% Tighten up the tolerances
sys.tols.rlx_integration=1e-5;

% Set relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={100e-12};

% Relaxation superoperator accuracy
sys.tols.rlx_integration=1e-5;
sys.tols.rlx_zero=1e-5;

% Use complete basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build the relaxation superoperator
R=relaxation(spin_system);

% List twenty smallest magnitude eigenvalues
disp('Twenty smallest relaxation rates, Hz:');
disp(eigs(R-speye(size(R)),20,'SM')+1);

end

