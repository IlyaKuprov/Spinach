% Relaxation analysis for strychnine, dipolar processes only.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function t1t2_strychnine()

% Spin system properties
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=5.9;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='kite';
inter.tau_c={200e-12};

% Distance cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation analysis
relaxan(spin_system);

end

