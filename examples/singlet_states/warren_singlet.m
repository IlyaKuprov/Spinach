% A demonstration that long-lived states exist that are immune
% not only to dipolar and CSA, but also to quadrupolar relaxati-
% on in certain circumstances. Full Redfield superoperator for
% dipolar and quadrupolar relaxation in liquid state is compu-
% ted and diagonalized.
%
% Calculation time: seconds
%
% wwarren@duke.edu
% i.kuprov@soton.ac.uk

function warren_singlet()

% System specification
sys.magnet=14.1;
sys.isotopes={'14N','14N'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.6 0.8 1.0]};
inter.coupling.matrix{1,1}=eeqq2nqi(1.25e6,0.25,1,[0 0 0]);
inter.coupling.matrix{2,2}=eeqq2nqi(1.25e6,0.25,1,[0 0 0]);

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={5e-9};

% Relaxation superoperator accuracy
sys.tols.rlx_integration=1e-5;
sys.tols.rlx_zero=1e-5;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);
sort(eig(full(R)))

end

