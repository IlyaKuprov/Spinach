% Tests operator expansion and conversion helpers. Syntax:
%
%                    result=test_operator_expansion_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that irreducible spherical tensor and bosonic monomial
% expansion helpers reconstruct explicit matrices, and that operator-sized
% allocation helpers return the correct formalism dimensions.
%
% ilya.kuprov@weizmann.ac.il

function result=test_operator_expansion_suite()

% State the expansion target of the test
result=new_test_result('kernel/operator_expansion_suite',...
                       'Operator expansions and conversions',...
                       'operator expansion coefficients must reconstruct the source matrices.');

% Check Hilbert-to-Liouville vectorisation identities on a non-diagonal matrix
H=[1 2;3 4];
unit=speye(2);
result=test_close(result,'hilb2liouv left explicit',hilb2liouv(H,'left'),kron(unit,H),1e-15,1e-15,...
                  'left multiplication vectorises as kron(I,H)');
result=test_close(result,'hilb2liouv right explicit',hilb2liouv(H,'right'),kron(transpose(H),unit),1e-15,1e-15,...
                  'right multiplication vectorises as kron(transpose(H),I)');
result=test_close(result,'hilb2liouv state vector explicit',hilb2liouv(H,'statevec'),H(:),1e-15,1e-15,...
                  'state-vector conversion stacks matrix columns');

% Check IST expansion of a generic spin-one matrix
A=[1 2 0;0 -1 3;4 0 2];
[states,coeffs]=oper2ist(A);
result=test_close(result,'oper2ist reconstruction',ist_reconstruct(3,states,coeffs),A,1e-12,1e-12,...
                  'complete rank-0, rank-1, and rank-2 tensors span all 3x3 spin-one matrices');

% Check spin and boson energy-level counting conventions in IST expansions
[states,coeffs]=enlev2ist(3,1,'S');
P=zeros(3); P(3,3)=1;
result=test_close(result,'enlev2ist spin bottom level',ist_reconstruct(3,states,coeffs),P,1e-12,1e-12,...
                  'spin energy levels are counted from the bottom upward');
[states,coeffs]=enlev2ist(3,1,'B');
P=zeros(3); P(1,1)=1;
result=test_close(result,'enlev2ist boson top level',ist_reconstruct(3,states,coeffs),P,1e-12,1e-12,...
                  'boson energy levels are counted from the top downward in IST conversion');

% Check central-transition and boson-product IST expansion wrappers
[states,coeffs]=ct2ist(4,'+');
result=test_close(result,'ct2ist reconstruction',ist_reconstruct(4,states,coeffs),centrans(4,'+'),1e-12,1e-12,...
                  'ct2ist() must expand the same central-transition matrix as centran() builds');
W=weyl(3);
[states,coeffs]=bos2ist('CAN',3);
result=test_close(result,'bos2ist reconstruction',ist_reconstruct(3,states,coeffs),W.c*W.a*W.n,1e-12,1e-12,...
                  'bos2ist() multiplies creation, annihilation, and number operators before IST expansion');

% Check bosonic monomial expansion of a generic finite oscillator matrix
A=[1 2 0;0 -1 3;4 0 2];
[states,coeffs]=oper2bm(A);
result=test_close(result,'oper2bm reconstruction',bm_reconstruct(3,states,coeffs),A,1e-11,1e-11,...
                  'finite oscillator monomials span all matrices in the truncated mode space');
[states,coeffs]=enlev2bm(3,2);
P=zeros(3); P(2,2)=1;
result=test_close(result,'enlev2bm reconstruction',bm_reconstruct(3,states,coeffs),P,1e-11,1e-11,...
                  'enlev2bm() expands the requested bosonic population projector');

% Build one-spin systems in the main formalisms
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.approximation='none';

% Check Hilbert-space unit and sparse preallocation dimensions
bas.formalism='zeeman-hilb';
spin_system=test_spin_system(sys,inter,bas);
result=test_close(result,'unit_oper zeeman-hilb',unit_oper(spin_system),speye(2),1e-15,1e-15,...
                  'Hilbert-space unit dimension is the spin multiplicity');
A=mprealloc(spin_system,2);
result=test_true(result,'mprealloc zeeman-hilb size',isequal(size(A),[2 2])&&(nnz(A)==0),...
                 'Hilbert-space preallocation creates an empty 2x2 sparse matrix');

% Check Zeeman-Liouville unit and sparse preallocation dimensions
bas.formalism='zeeman-liouv';
spin_system=test_spin_system(sys,inter,bas);
result=test_close(result,'unit_oper zeeman-liouv',unit_oper(spin_system),speye(4),1e-15,1e-15,...
                  'Zeeman-Liouville unit dimension is the square of the spin multiplicity');
A=mprealloc(spin_system,2);
result=test_true(result,'mprealloc zeeman-liouv size',isequal(size(A),[4 4])&&(nnz(A)==0),...
                 'Zeeman-Liouville preallocation creates an empty 4x4 sparse matrix');

% Check spherical-tensor Liouville dimensions without assuming a fixed basis size literal
bas.formalism='sphten-liouv';
spin_system=test_spin_system(sys,inter,bas);
basis_dim=size(spin_system.bas.basis,1);
result=test_close(result,'unit_oper sphten-liouv',unit_oper(spin_system),speye(basis_dim),1e-15,1e-15,...
                  'spherical-tensor Liouville unit dimension is the number of retained basis states');
A=mprealloc(spin_system,2);
result=test_true(result,'mprealloc sphten-liouv size',isequal(size(A),[basis_dim basis_dim])&&(nnz(A)==0),...
                 'spherical-tensor preallocation follows the current basis size');

% Check two-spin irreducible tensor formula in Hilbert space
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{2,2}=0;
bas.formalism='zeeman-hilb';
spin_system=test_spin_system(sys,inter,bas);
T20=twospinist(spin_system,1,2,[2 0],'comm');
T20_ref=sqrt(2/3)*(operator(spin_system,{'Lz','Lz'},{1,2})-...
        (1/4)*(operator(spin_system,{'L+','L-'},{1,2})+...
               operator(spin_system,{'L-','L+'},{1,2})));
result=test_close(result,'twospinist rank two zero projection',T20,T20_ref,1e-15,1e-15,...
                  'the two-spin rank-two zero-projection tensor follows the Clebsch-Gordan combination');

% Check Lindbladian rate calibration on a diagonal jump process
rho=[1;1];
rate=2.75;
R=lindbladian(diag([1 0]),diag([0 1]),rho,rate);
obs=real((rho'*R*rho)/(rho'*rho));
result=test_close(result,'lindbladian requested rate',obs,-rate,1e-12,1e-12,...
                  'lindbladian() rescales the generator to the requested decay rate of rho');

end


function A=ist_reconstruct(mult,states,coeffs)

% Build the complete IST basis for the requested multiplicity
T=irr_sph_ten(mult);

% Recombine zero-based Spinach IST indices into a matrix
A=zeros(mult);
for n=1:numel(states)
    A=A+coeffs(n)*T{states(n)+1};
end

end


function A=bm_reconstruct(nlevels,states,coeffs)

% Build the complete bosonic monomial basis for the requested truncation
B=boson_mono(nlevels);

% Recombine zero-based Spinach BM indices into a matrix
A=zeros(nlevels);
for n=1:numel(states)
    A=A+coeffs(n)*B{states(n)+1};
end

end


