% A demonstration that the nitrogen singlet state in urea is not
% long-lived. The relaxation superoperator accounts for every di-
% polar coupling and every CSA tensor in the system.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function decoherence_urea()

% Read the spin system (coordinates, chemical shifts,
% J-couplings and CSAs) from a vacuum DFT calculation
[sys,inter]=g2spinach(gparse('../standard_systems/urea.log'),...
                   {{'H','1H'},{'N','15N'}},[30.0 166.0],[]);

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

% Get the Lz and the singlet for 15N
Lz=state(spin_system,'Lz','15N');
S=singlet(spin_system,1,4);

% See what remains of them under R
disp(['Norm of R|Lz>: ' num2str(norm(R*Lz)/norm(Lz))]);
disp(['Norm of R|S>:  ' num2str(norm(R*S)/norm(S))]);

end

