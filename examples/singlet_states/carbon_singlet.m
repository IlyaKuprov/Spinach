% Singlet relaxation rate for the two triple bond carbons
% in cis-dimethylbut-2-ynedioate. Magnetic parameters com-
% puted with DFT. 
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function carbon_singlet()

% System specification
sys.magnet=14.1;
sys.isotopes={'13C','13C'};
inter.zeeman.matrix{1}=[29.13        0.00        0.00
                         0.00      253.00       -5.64
                         0.00      -35.61       38.78];
inter.zeeman.matrix{2}=[29.13        0.00        0.00
                         0.00      253.00        5.64
                         0.00       35.61       38.78];
inter.coordinates={[0.000    0.609    0.298]
                   [0.000   -0.609    0.298]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={100e-12};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation superoperator accuracy
sys.tols.rlx_integration=1e-5;
sys.tols.rlx_zero=1e-5;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Action on longitudinal magnetization
Sz=state(spin_system,{'Lz'},{1}); Sz=Sz/norm(Sz);
report(spin_system,['<Sz|R|Sz> matrix element: ' num2str(Sz'*R*Sz)]);

% Action on a singlet state
S=singlet(spin_system,1,2); S=S/norm(S);
report(spin_system,['<singlet|R|singlet> matrix element: ' num2str(S'*R*S)]);

end

