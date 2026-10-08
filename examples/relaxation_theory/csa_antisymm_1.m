% Longitudinal and transverse relaxation rates in a system 
% with a significant antisymmetry in the shielding tensor.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function csa_antisymm_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'13C'};
inter.zeeman.matrix={[100  20  15
                      20   0   30
                      25   10 -30]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={50e-12};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Spinach relaxation rates 
R=relaxation(spin_system);
Lz=state(spin_system,{'Lz'},{1});
Lp=state(spin_system,{'L+'},{1});
R1Sp=-(Lz'*R*Lz)/(Lz'*Lz);
R2Sp=-(Lp'*R*Lp)/(Lp'*Lp);

% Textbook relaxation rates
[R1Book,R2Book]=rlx_csa(sys.magnet,sys.isotopes{1},...
                inter.zeeman.matrix{1},inter.tau_c{1});

% Summary
disp(['R1 rate, Spinach:  ' num2str(R1Sp)]);
disp(['R1 rate, textbook: ' num2str(R1Book)]);
disp(['R2 rate, Spinach:  ' num2str(R2Sp)]);
disp(['R2 rate, textbook: ' num2str(R2Book)]);
   
end

