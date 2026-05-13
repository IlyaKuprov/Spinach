% Tests Hilbert-space operator generation. Syntax:
%
%                    result=test_hilbert_operator()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that operator() builds the correct one-spin Hilbert-space
% angular momentum matrices from human-readable labels.
%
% ilya.kuprov@weizmann.ac.il

function result=test_hilbert_operator()

% Announce the test target
fprintf('TESTING: Hilbert-space operator generation\n');

% State the physical target of the test
result=new_test_result('kernel/hilbert_operator',...
                       'Hilbert-space operator generation',...
                       'operator() must map Lx, Ly, Lz, L+, and L- labels to spin matrices.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Textbook spin-half reference matrices
S=pauli(2);

% Check label-to-matrix mapping
result=test_close(result,'Lx operator',operator(spin_system,'Lx',1),S.x,1e-15,1e-15,...
                  'a one-proton Lx operator is the spin-half Sx matrix');
result=test_close(result,'Ly operator',operator(spin_system,'Ly',1),S.y,1e-15,1e-15,...
                  'a one-proton Ly operator is the spin-half Sy matrix');
result=test_close(result,'Lz operator',operator(spin_system,'Lz',1),S.z,1e-15,1e-15,...
                  'a one-proton Lz operator is the spin-half Sz matrix');
result=test_close(result,'L+ operator',operator(spin_system,'L+',1),S.p,1e-15,1e-15,...
                  'L+ is the spin raising operator in the Zeeman basis');
result=test_close(result,'L- operator',operator(spin_system,'L-',1),S.m,1e-15,1e-15,...
                  'L- is the spin lowering operator in the Zeeman basis');

end

