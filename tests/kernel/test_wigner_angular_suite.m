% Tests angular-momentum coefficient and spherical-function helpers. Syntax:
%
%                    result=test_wigner_angular_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Clebsch-Gordan coefficients, Wigner symbols, Wigner D
% matrices, and spherical harmonics against elementary exact values.
%
% ilya.kuprov@weizmann.ac.il

function result=test_wigner_angular_suite()

% Announce the test target
fprintf('TESTING: Angular-momentum coefficient functions\n');

% State the angular target of the test
result=new_test_result('kernel/wigner_angular_suite',...
                       'Angular-momentum coefficient functions',...
                       'angular-momentum helpers must reproduce elementary exact quantum-mechanical coefficients.');

% Coupling two spin-half particles gives triplet and singlet M=0 amplitudes of 1/sqrt(2)
result=test_close(result,'clebsch_gordan triplet half-half',...
                  clebsch_gordan(1,0,1/2,1/2,1/2,-1/2),1/sqrt(2),1e-14,1e-14,...
                  'the |1,0> triplet contains |alpha beta> with coefficient 1/sqrt(2)');
result=test_close(result,'clebsch_gordan singlet half-half',...
                  clebsch_gordan(0,0,1/2,1/2,1/2,-1/2),1/sqrt(2),1e-14,1e-14,...
                  'the |0,0> singlet contains |alpha beta> with coefficient 1/sqrt(2) in Spinach phase convention');
result=test_close(result,'clebsch_gordan forbidden projection',...
                  clebsch_gordan(1,1,1/2,1/2,1/2,-1/2),0,1e-14,1e-14,...
                  'forbidden projection combinations return zero');

% Wigner 3j values follow from the relation to Clebsch-Gordan coefficients
result=test_close(result,'wigner_3j 1 1 0',wigner_3j(1,0,1,0,0,0),-1/sqrt(3),1e-14,1e-14,...
                  'the elementary (1 1 0; 0 0 0) Wigner 3j symbol is -1/sqrt(3)');
result=test_close(result,'wigner_6j all zero',wigner_6j(0,0,0,0,0,0),1,1e-14,1e-14,...
                  'the all-zero Wigner 6j symbol is unity');
result=test_close(result,'wigner_6j all ones',wigner_6j(1,1,1,1,1,1),1/6,1e-14,1e-14,...
                  'the Wigner 6j symbol with all angular momenta one is 1/6');

% Wigner D matrices are unitary representations and reduce to identity for zero rotation
D1=wigner(1,0,pi/2,0);
D1_ref=[1/2 -1/sqrt(2) 1/2; 1/sqrt(2) 0 -1/sqrt(2); 1/2 1/sqrt(2) 1/2];
result=test_close(result,'wigner rank one beta pi/2',D1,D1_ref,1e-14,1e-14,...
                  'the rank-one Wigner matrix at beta=pi/2 has the Brink-Satchler closed form');
D0=wigner(2,0,0,0);
D=wigner(2,0.2,0.4,0.7);
result=test_close(result,'wigner identity rotation',D0,eye(5),1e-14,1e-14,...
                  'zero Euler angles give the identity Wigner D matrix');
result=test_close(result,'wigner unitarity',D'*D,eye(5),1e-13,1e-13,...
                  'Wigner D matrices are unitary rotation representations');

% Spherical harmonics have elementary normalised values
th=[0 pi/2 pi]; ph=[0 pi/3 pi/7];
result=test_close(result,'spher_harmon Y00',spher_harmon(0,0,th,ph),ones(size(th))/sqrt(4*pi),1e-14,1e-14,...
                  'Y_0^0 is the constant 1/sqrt(4*pi)');
result=test_close(result,'spher_harmon Y10',spher_harmon(1,0,th,ph),sqrt(3/(4*pi))*cos(th),1e-14,1e-14,...
                  'Y_1^0 is sqrt(3/(4*pi))*cos(theta)');

end
