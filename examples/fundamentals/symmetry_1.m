% Liouvillian symmetrization for a radical pair with four 
% equivalent nuclei under the S4 permutation group.
%
% ilya.kuprov@weizmann.ac.il

function symmetry_1()

% Magnetic field
sys.magnet=0;

% Spin system
sys.isotopes ={'E','E','1H','1H','1H','1H'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_spins={[3 4 5 6]};
bas.sym_group={'S4'};
bas.sym_a1g_only=0;

% Interactions
inter.zeeman.scalar={2.002 2.002 0 0 0 0};
inter.coupling.scalar=num2cell(mt2hz([0      0   0.295 0.295 0.295 0.295
                                      0      0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0]));
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'labframe');

% Hamiltonian superoperator
H=hamiltonian(spin_system);

% Symmetry factorization
S=horzcat(spin_system.bas.irrep.projector);

% Plotting
kfigure(); scale_figure([1.5 1]);
subplot(1,2,1); spy(abs(H)>1e3); 
ktitle('Original Liouvillian');
subplot(1,2,2); spy(abs(S'*H*S)>1e3); 
ktitle('Symmetrised Liouvillian');
xline(560); yline(560);
xline(1216); yline(1216); xline(1936); yline(1936);

end

