% A demonstration that the two-spin singet state is immune 
% to dipolar relaxation. Full Redfield superoperator for
% dipolar relaxation in liquid state is computed and the
% norm of its action on a singlet state is printed to the
% console.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dipolar_singlet()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.5 0.6 0.7]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={5e-9};

% Relaxation superoperator accuracy
sys.tols.rlx_integration=1e-5;
sys.tols.rlx_zero=1e-5;

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Action on a singlet state
S=singlet(spin_system,1,2); S=S/norm(S);
disp(['Norm of R|S>: ' num2str(norm(R*S))]);

end

