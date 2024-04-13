% Complete Bloch-Redfield-Wangsness relaxation superoperator in a system 
% with a dipolar coupling between two spins. Dipolar couplings are compu-
% ted from Cartesian coordinates of the two spins.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function dd_relaxation_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 1.02]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={5e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Textbook rates
[r1,r2,rx]=rlx_dip(sys.magnet,sys.isotopes,...
                   norm(inter.coordinates{1}-...
                        inter.coordinates{2}),inter.tau_c{1});

% Textbook and Spinach R1 for first spin
rho=state(spin_system,'Lz',1); rho=rho/norm(rho,2);
disp(['(Spin 1, ' sys.isotopes{1} ') R1 textbook: '       ...
                  num2str(r1(1),'%10.6f') ' Hz, Spinach: ' ...
                  num2str(-rho'*R*rho,'%10.6f') ' Hz']);

% Textbook and Spinach R1 for second spin
rho=state(spin_system,'Lz',2); rho=rho/norm(rho,2);
disp(['(Spin 2, ' sys.isotopes{2} ') R1 textbook: '       ...
                  num2str(r1(2),'%10.6f') ' Hz, Spinach: ' ...
                  num2str(-rho'*R*rho,'%10.6f') ' Hz']);

% Textbook and Spinach R2 for first spin
rho=state(spin_system,'L+',1); rho=rho/norm(rho,2);
disp(['(Spin 1, ' sys.isotopes{1} ') R2 textbook: '       ...
                  num2str(r2(1),'%10.6f') ' Hz, Spinach: ' ...
                  num2str(-rho'*R*rho,'%10.6f') ' Hz']);

% Textbook and Spinach R2 for second spin
rho=state(spin_system,'L+',2); rho=rho/norm(rho,2);
disp(['(Spin 2, ' sys.isotopes{2} ') R2 textbook: '       ...
                  num2str(r2(2),'%10.6f') ' Hz, Spinach: ' ...
                  num2str(-rho'*R*rho,'%10.6f') ' Hz']);

% Print textbook cross-relaxation rate
rho_a=state(spin_system,'Lz',1); rho_a=rho_a/norm(rho_a,2);
rho_b=state(spin_system,'Lz',2); rho_b=rho_b/norm(rho_b,2);
disp(['(1) -> (2) Rx, textbook: ' num2str(rx,'%10.6f') ' Hz, Spinach: ' ...
                                  num2str(-rho_b'*R*rho_a,'%10.6f') ' Hz']);

% Print the full relaxation superoperator
disp('Complete relaxation superoperator, IST basis:'); disp(full(R));

end

