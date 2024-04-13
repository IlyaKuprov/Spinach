% Test of matrix exponential differentiation routines. The first pair of
% matrices must be equal; the second pair of matrices must be equal.
%
% i.kuprov@soton.ac.uk

function dirdiff_1()

% Isotopes
sys.isotopes={'1H','1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={1.0 1.5};

% Scalar couplings
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=7.0; 

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Random Hamiltonian
H=randn(5)+1i*randn(5); H=(H+H')/2;

% Random direction operators
A=randn(5)+1i*randn(5); A=(A+A')/20;
B=randn(5)+1i*randn(5); B=(B+B')/20;

% First derivative, numerical
D=(propagator(spin_system,H+1e-3*A,1)-propagator(spin_system,H-1e-3*A,1))/2e-3; 
disp('First derivative, numerical:'); disp(D); disp(' ');

% First derivative, analytical
D=dirdiff(spin_system,H,A,1,2);
disp('First derivative, analytical:'); disp(D{2}); disp(' ');

% Second derivative, numerical
D=(propagator(spin_system,H+1e-3*A+1e-3*B,1)-propagator(spin_system,H+1e-3*A-1e-3*B,1)-...
   propagator(spin_system,H-1e-3*A+1e-3*B,1)+propagator(spin_system,H-1e-3*A-1e-3*B,1))/4e-6;
disp('Second derivative, numerical:'); disp(D); disp(' ');

% Second derivative, analytical
P=dirdiff(spin_system,H,{A,B},1,3);
Q=dirdiff(spin_system,H,{B,A},1,3);
disp('Second derivative, analytical:');
disp((P{3}+Q{3})/2);

end

