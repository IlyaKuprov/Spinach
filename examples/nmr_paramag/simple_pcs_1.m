% Pseudocontact shift and Curie relaxation on a proton due to
% the presence of a point magnetic susceptibility centre.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function simple_pcs_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H'};

% Diamagnetic shifts and coordinates
inter.zeeman.scalar={2.0};
inter.coordinates={[0.0 0.0 0.0]};

% Magnetic susceptibility
inter.suscept.chi={[0.0883   -0.0904    0.0822
                   -0.0904   -0.1011   -0.0149
                    0.0822   -0.0149    0.0128]};
inter.suscept.xyz={[10.0 2.5 3.9]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={10e-12};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Spinach relaxation rates 
R=relaxation(spin_system);
Lz=state(spin_system,'Lz','1H');
Lp=state(spin_system,'L+','1H');
R1Sp=-(Lz'*R*Lz)/(Lz'*Lz);
R2Sp=-(Lp'*R*Lp)/(Lp'*Lp);
      
% Summary
disp(['R1 rate, Spinach: ' num2str(R1Sp)]);
disp(['R2 rate, Spinach: ' num2str(R2Sp)]);

end

