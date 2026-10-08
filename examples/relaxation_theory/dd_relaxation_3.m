% Extreme narrowing limit case comparison between the dipolar relaxation
% rates in proton-proton and proton-deuterium system. The rate must sca-
% le with the square of the magnetogyric ratio and with S(S+1), where S
% is the quantum number of the partner spin.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dd_relaxation_3()

%% System with two protons

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 2.50]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1e-12};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Longitudinal relaxation rate
rho=state(spin_system,{'Lz'},{1}); 
rho=rho/norm(rho); R1pp=rho'*R*rho;

%% System with a proton and a deuteron

% System specification - two protons
sys.magnet=14.1;
sys.isotopes={'1H','2H'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 2.50]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1e-12};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Longitudinal relaxation rate
rho=state(spin_system,{'Lz'},{1}); 
rho=rho/norm(rho); R1pd=rho'*R*rho;

%% Comparison
disp(['R1 rate ratio, Spinach: ' num2str(R1pp/R1pd)]);
ratio=(((1/2)*(1/2+1))/(1*(1+1)))*(spin('1H')/spin('2H'))^2;
disp(['R1 rate ratio, textbook: ' num2str(ratio)]);

end

