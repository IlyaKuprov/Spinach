% Test of the reverse decomposition of spin-1 
% Hamiltonians, useful for transmon parameter
% import.
%
% ilya.kuprov@weizmann.ac.il

function transmon_import()

% Get a random three-level
% transmon Hamiltonian
H_T=randn(3,3)+1i*randn(3,3); 
H_T=H_T+H_T'; H_T=remtrace(H_T);

% Translate into Spinach
[omega,Q]=ham2nqi(H_T);

% Set up Spinach
sys.magnet=0;
sys.isotopes={'T3'};
inter.coupling.matrix={Q/(2*pi)};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Re-build using Spinach functionality
spin_system=assume(spin_system,'labframe');
[H_iso,H_aniso]=hamiltonian(spin_system);
H_S=H_iso+orientation(H_aniso,[0 0 0]);
H_S=H_S+omega(1)*operator(spin_system,'Lx','T3');
H_S=H_S+omega(2)*operator(spin_system,'Ly','T3');
H_S=H_S+omega(3)*operator(spin_system,'Lz','T3');

% Compare the matrices
disp(H_T); disp(' '); disp(full(H_S));

end