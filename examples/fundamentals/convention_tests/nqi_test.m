% Test of the reverse decomposition of
% spin-1 Hamiltonians.
%
% ilya.kuprov@weizmann.ac.il

function nqi_test()

% Get random test Hamiltonian
H_T=randn(3,3)+1i*randn(3,3); 
H_T=H_T+H_T'; H_T=remtrace(H_T);

% Translate back
[omega,Q]=ham2nqi(H_T);

% Set up Spinach
sys.magnet=0;
sys.isotopes={'14N'};
inter.coupling.matrix={Q/(2*pi)};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Re-build using Spinach functionality
spin_system=assume(spin_system,'labframe');
[H_iso,H_aniso]=hamiltonian(spin_system);
H_S=H_iso+orientation(H_aniso,[0 0 0]);
H_S=H_S+omega(1)*operator(spin_system,'Lx','14N');
H_S=H_S+omega(2)*operator(spin_system,'Ly','14N');
H_S=H_S+omega(3)*operator(spin_system,'Lz','14N');

% Compare the matrices
if norm(H_T-H_S,2)>1e-6*norm(H_T+H_S,2)
    error('Quadrupolar reconstruction test FAILED.');
else
    disp('Quadrupolar reconstruction test PASSED.');
end

end