% Tests spin-one angular momentum matrices. Syntax:
%
%                    result=test_pauli_spin_one_algebra()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the spin-one representation: Sz projections are +1, 0,
% and -1, ladder matrix elements are sqrt(2), and S^2=s(s+1)=2.
%
% ilya.kuprov@weizmann.ac.il

function result=test_pauli_spin_one_algebra()

% State the physical target of the test
result=new_test_result('kernel/pauli_spin_one_algebra',...
                       'Spin-one angular momentum algebra',...
                       'Spin-one operators must realise su(2), with S^2=2.');

% Generate Spinach spin-one operators
S=pauli(3);

% Write the textbook spin-one matrices explicitly
Sp=[0 sqrt(2) 0;0 0 sqrt(2);0 0 0];
Sm=Sp';
Sz=diag([1 0 -1]);
Sx=(Sp+Sm)/2;
Sy=(Sp-Sm)/2i;

% Check matrix elements and commutators
result=test_close(result,'Sz projections',S.z,Sz,1e-15,1e-15,...
                  'spin-one magnetic quantum numbers are +1, 0, and -1');
result=test_close(result,'S+ matrix elements',S.p,Sp,1e-15,1e-15,...
                  'ladder matrix elements are sqrt(s(s+1)-m(m+1))');
result=test_close(result,'Sx from ladder operators',S.x,Sx,1e-15,1e-15,...
                  'Cartesian Sx is the Hermitian half-sum of S+ and S-');
result=test_close(result,'Sy from ladder operators',S.y,Sy,1e-15,1e-15,...
                  'Cartesian Sy is the Hermitian half-difference divided by i');
result=test_close(result,'[Sx,Sy]=iSz',comm(S.x,S.y),1i*S.z,1e-15,1e-15,...
                  'the spin-one matrices must obey the same su(2) algebra');

% Check the Casimir operator
S2=S.x*S.x+S.y*S.y+S.z*S.z;
result=test_close(result,'S^2=s(s+1)',S2,2*S.u,1e-15,1e-15,...
                  'for s=1 the Casimir eigenvalue is 2');

end

