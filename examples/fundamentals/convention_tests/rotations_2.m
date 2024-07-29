% A rotations test comparing the Hamiltonians for a manually rotated (at
% the interaction specification level) spin system with the Hamiltonian
% that has been rotated using Spinach operator rotation functionality.
%
% i.kuprov@soton.ac.uk
% john.price@colorado.edu

function rotations_2()

% Generate random matrices
shift_tensor_a=randn(3,3);
shift_tensor_b=randn(3,3);
coupl_tensor_a=100*randn(3,3);

%% Kernel level rotation

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% A pair of spins at a distance, A
sys.isotopes={'1H','15N'};
inter.zeeman.matrix={shift_tensor_a...
                     shift_tensor_b};
inter.coordinates={[0.7 0.8 0.9];
                   [1.5 2.5 3.5]};
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=coupl_tensor_a;
               
% Spinach housekeeping, A
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian, A
[H,Q]=hamiltonian(assume(spin_system,'labframe'));
H_A=H+orientation(Q,[1.0 2.0 3.0]);

%% Input level rotation

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% A pair of spins at a distance, B
sys.isotopes={'1H','15N'};
R=euler2dcm(1.0,2.0,3.0);
inter.zeeman.matrix={R*shift_tensor_a*R'...
                     R*shift_tensor_b*R'};
inter.coordinates={[0.7 0.8 0.9]*R';
                   [1.5 2.5 3.5]*R'};
inter.coupling.matrix{1,2}=R*coupl_tensor_a*R';

% Spinach housekeeping, B
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian, B
[H,Q]=hamiltonian(assume(spin_system,'labframe'));
H_B=H+orientation(Q,[0.0 0.0 0.0]);

%% Comparison

% Norm of the residual
resnorm=norm(H_A-H_B,1);

% Display the diagnostics
if resnorm>1e-3
    report(spin_system,['rotation test failed, ' num2str(resnorm)]);
else
    report(spin_system,['rotation test passed, ' num2str(resnorm)]);
end

end

