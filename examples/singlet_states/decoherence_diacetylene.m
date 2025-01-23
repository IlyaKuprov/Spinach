% Long-lived spin states in the diacetylene molecule (2 protons, 
% 4 carbons, 4096-dimensional Liouville space). The relaxation
% superoperator accounts for every dipolar coupling and every CSA
% tensor in the system.
%
% Calculation time: seconds
%
% yesu.feng@duke.edu
% ilya.kuprov@weizmann.ac.il

function decoherence_diacetylene()

% Read the spin system (coordinates, chemical shifts,
% J-couplings and CSAs) from a vacuum DFT calculation
[sys,inter]=g2spinach(gparse('../standard_systems/diacetylene.log'),...
                          {{'H','1H'},{'C','13C'}},[31.8 182.4],[]);

% Set magnet field to 1.0 Tesla
sys.magnet=14.1;

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
disp(' '); format('long');

% List twenty smallest magnitude eigenvalues of the relaxation
% superoperator (diagonal preconditioning is used)
disp('20 smallest eigenvalues of the relaxation superoperator:');
disp(eigs(R-speye(size(R)),20,'SM')+1);

% Compute the self-relaxation rate of the singlet state between the two
% centre carbons (spins 1 and 2 in this case)
S=singlet(spin_system,1,2); S=S/norm(S);
disp('Self-relaxation rate of the two-carbon singlet:');
disp(S'*R*S);

% Find the eigenvectors corresponding to the slowly relaxing states
[v,~]=eigs(R-speye(size(R)),2,'SM');

% Get the spherical tensor composition of the slowly relaxing states
disp('Slowly relaxing state 1:');
stateinfo(spin_system,v(:,1),100);
disp('Slowly relaxing state 2:');
stateinfo(spin_system,v(:,2),100);

end

