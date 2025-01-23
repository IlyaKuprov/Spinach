% Complete Bloch-Redfield-Wangsness relaxation superoperator in a system 
% with two anisotropically shielded nuclei. Spinach relaxation theory mo-
% dule automatically accounts for all cross-correlations (CSA-CSA cross-
% correlation is present in this case).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function csa_csa_xcorr_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.eigs={[7 15 -22]
                   [11 18 -29]};
inter.zeeman.euler={[pi/5 pi/3 pi/11]
                    [pi/6 pi/7 pi/15]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={2e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
disp(full(relaxation(spin_system)));

end

