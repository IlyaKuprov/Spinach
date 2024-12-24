% Complete Bloch-Redfield-Wangsness relaxation superoperator in a system 
% with two anisotropically shielded nuclei with a dipolar coupling betwe-
% en them. Spinach relaxation theory module automatically accounts for all
% cross-correlations (CSA-CSA and DD-CSA cross-correlations are both pre-
% sent in this case). Dipolar couplings are computed from Cartesian coor-
% dinates of the two spins.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dd_csa_xcorr_1()

% Spin system
sys.isotopes={'1H','13C'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Interactions
sys.magnet=14.1;
inter.zeeman.eigs={[7  15 -22]
                   [11 18 -29]};
inter.zeeman.euler={[pi/3 pi/4 pi/5]
                    [pi/6 pi/7 pi/8]};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.02]};
               
% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1e-9};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
disp(full(relaxation(spin_system)));

end

