% Complete Bloch-Redfield-Wangsness relaxation superoperator in a system 
% with dipolar coupling between spins. The dipolar couplings are computed
% from Cartesian coordinates of the spins. The result should not depend
% on the choice of the rotation angles below.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dd_relaxation_2()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={1.0,2.0,3.0};

% Randomly rotated set of coordinates
R=euler2dcm([pi/3 pi/4 pi/5]);
inter.coordinates={[-0.2230    0.7893   -0.5721]*R
                   [-0.6752   -0.3561    0.6460]*R
                   [ 0.8982   -0.4333   -0.0739]*R};

% Relaxation theory parameters
sys.tols.rlx_integration=1e-5;
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
disp(full(relaxation(spin_system)));

end

