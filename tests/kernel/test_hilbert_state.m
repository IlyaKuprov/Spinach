% Tests Hilbert-space state generation. Syntax:
%
%                    result=test_hilbert_state()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that state() returns the expected density matrices in
% Hilbert space for a one-spin system.
%
% ilya.kuprov@weizmann.ac.il

function result=test_hilbert_state()

% State the physical target of the test
result=new_test_result('kernel/hilbert_state',...
                       'Hilbert-space state generation',...
                       'state() must map observable labels to density matrices.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Textbook spin-half reference matrices
S=pauli(2);

% Check density matrices generated from state labels
result=test_close(result,'Lz state',state(spin_system,'Lz',1),S.z,1e-15,1e-15,...
                  'the Lz state is the longitudinal magnetisation density matrix');
result=test_close(result,'Lx state',state(spin_system,'Lx',1),S.x,1e-15,1e-15,...
                  'the Lx state is the transverse in-phase density matrix');
result=test_close(result,'identity state',state(spin_system,'E',1),S.u,1e-15,1e-15,...
                  'the E state is the unit density matrix');

end

