% Deuterium pair singlet, triplet, and quintet state
% internal consistency test.
%
% ilya.kuprov@weizmann.ac.il

function state_consistency_3()

% A pair of deuteria
sys.magnet=0;
sys.isotopes={'2H','2H'};
inter.zeeman.scalar={0 0};

% Hilbert space
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Ortho-deuterium states from Spinach
[S,T,Q]=deut_pair(spin_system,1,2);

% Component vectors in Hilbert space as per Eq 1
% in https://doi.org/10.1016/S0009-2614(98)00784-2
alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1];
S0=(1/sqrt(3))*(kron(alp,gam)-kron(bet,bet)+kron(gam,alp));
Tp=(1/sqrt(2))*(kron(alp,bet)-kron(bet,alp));
T0=(1/sqrt(2))*(kron(alp,gam)-kron(gam,alp));
Tm=(1/sqrt(2))*(kron(bet,gam)-kron(gam,bet));
Qpp=kron(alp,alp); 
Qp=(1/sqrt(2))*(kron(alp,bet)+kron(bet,alp));
Q0=(1/sqrt(6))*(kron(alp,gam)+2*kron(bet,bet)+kron(gam,alp));
Qm=(1/sqrt(2))*(kron(bet,gam)+kron(gam,bet));
Qmm=kron(gam,gam);

% Test singlet state
if norm(S-S0*S0',1)>1e-6
    error('State construction test FAILED.');
end

% Test triplet states
if (norm(T{1}-Tp*Tp',1)>1e-6)||...
   (norm(T{2}-T0*T0',1)>1e-6)||...
   (norm(T{3}-Tm*Tm',1)>1e-6)
    error('State construction test FAILED.');
end

% Test quintet states
if (norm(Q{1}-Qpp*Qpp',1)>1e-6)||...
   (norm(Q{2}-Qp*Qp',1)>1e-6)||...
   (norm(Q{3}-Q0*Q0',1)>1e-6)||...
   (norm(Q{4}-Qm*Qm',1)>1e-6)||...
   (norm(Q{5}-Qmm*Qmm',1)>1e-6)
    error('State construction test FAILED.');
end

% Move to Liouville space
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Ortho-deuterium states from Spinach
[S,T,Q]=deut_pair(spin_system,1,2);

% Unit state from Spinach
U=state(spin_system,{'E','E'},{1 2});

% Summation into the unit state test
if norm(S+T{1}+T{2}+T{3}+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}-U,1)>1e-6
    error('State construction test FAILED.');
end

% Report success
disp('State construction test PASSED.');

end

