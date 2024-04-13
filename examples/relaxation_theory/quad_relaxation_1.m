% 14N quadrupolar relaxation in glycine in liquid state. The
% numerical output of Spinach is compared to the analytical
% equation from the textbook.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function quad_relaxation_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'14N'};

% Spin quantum number and quadrupolar tensor
[~,s_mult]=spin(sys.isotopes{1}); s_qnum=(s_mult-1)/2;
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,s_qnum,[0 0 0]);

% Relaxation theory
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
R=relaxation(spin_system);

% Textbook relaxation rate expressions
[r1,r2]=rlx_nqi(s_qnum,14.1*spin(sys.isotopes{1}),1.18e6,0.53,1e-9);

% States of interest
Lz=state(spin_system,'Lz',1);
Lp=state(spin_system,'L+',1);
Lz=Lz/norm(Lz,2); Lp=Lp/norm(Lp,2);

% Print the answers
disp([sys.isotopes{1} ' longitudinal relaxation rate, Spinach:  ' num2str(-Lz'*R*Lz)]);
disp([sys.isotopes{1} ' longitudinal relaxation rate, textbook: ' num2str(r1)]);
disp([sys.isotopes{1} ' transverse relaxation rate, Spinach:    ' num2str(-Lp'*R*Lp)]);
disp([sys.isotopes{1} ' transverse relaxation rate, textbook:   ' num2str(r2)]);

end

