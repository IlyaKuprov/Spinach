% Tests cheap deterministic low-level utility functions. Syntax:
%
%                    result=test_lowlevel_utilities_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks small numerical helpers, matrix filters, integer type
% selection, and analytic line-shape definitions against explicit answers.
%
% ilya.kuprov@weizmann.ac.il

function result=test_lowlevel_utilities_suite()

% State the utility target of the test
result=new_test_result('kernel/lowlevel_utilities_suite',...
                       'Low-level utility functions',...
                       'cheap scalar and matrix helpers must preserve their documented algebraic definitions.');

% Define small test matrices
A=[1 2;3 4];
B=[0 1;-1 2];

% Check commutator and right-ordered nested commutator
result=test_close(result,'comm',comm(A,B),A*B-B*A,1e-15,1e-15,...
                  'the commutator is AB-BA');
result=test_close(result,'rocomm',rocomm({A,B,A}),comm(comm(A,B),A),1e-15,1e-15,...
                  'rocomm nests commutators from left to right over the supplied cell array');

% Check trace removal and commuting part extraction
C=[2 1;0 4];
result=test_close(result,'remtrace',remtrace(C),[-1 1;0 1],1e-15,1e-15,...
                  'remtrace subtracts trace(A)/dim times the unit matrix');
H=[2 1+2i;1-2i 5];
result=test_close(result,'remncomm diagonal basis',remncomm(H,eye(2)),diag(diag(H)),1e-15,1e-15,...
                  'in the eigenbasis of B only the diagonal part commutes with B');

% Check Frobenius inner product and anti-diagonal transpose
D=[1+1i 2-1i;3 4i];
E=[2 0;1i -3];
result=test_close(result,'hdot',hdot(D,E),trace(D'*E),1e-15,1e-15,...
                  'hdot is the Frobenius inner product trace(A''B)');
result=test_close(result,'atranspose',atranspose(A),[4 2;3 1],1e-15,1e-15,...
                  'atranspose reflects a matrix across the anti-diagonal');

% Check matrix wiping helpers
M=reshape(1:16,4,4);
result=test_close(result,'killcross',killcross(M,[2 4],[1 3]),...
                  [0 0 0 0;2 0 10 0;0 0 0 0;4 0 12 0],1e-15,1e-15,...
                  'killcross zeroes the requested columns and rows');
result=test_close(result,'killdiag',killdiag(M,1),...
                  [0 0 9 13;0 0 0 14;3 0 0 0;4 8 0 0],1e-15,1e-15,...
                  'killdiag zeroes the requested diagonal brush band');

% Check rank and SVD truncation helpers
S=diag([5 2 1]);
result=test_close(result,'keep_rank',keep_rank(S,2),diag([5 2 0]),1e-14,1e-14,...
                  'keep_rank reconstructs a matrix from the requested leading singular values');
result=test_close(result,'frob_chop',frob_chop([5 4 0.3 0.1],0.5),2,1e-15,1e-15,...
                  'frob_chop keeps the smallest rank whose discarded tail is below tolerance');

% Check analytic line shapes and spectral density at points with closed-form values
x=[0 1];
g_fwhm=2*sqrt(2*log(2));
result=test_close(result,'gaussfun',gaussfun(x,g_fwhm),exp(-(x.^2)/2)/sqrt(2*pi),1e-15,1e-15,...
                  'with sigma one, gaussfun is the standard normal density');
[lor_r,lor_i]=lorentzfun(0,2*pi,2,x,0);
result=test_close(result,'lorentzfun real',lor_r,[1 1/2],1e-15,1e-15,...
                  'zero-phase Lorentzian real part is 1/(1+x^2) for the chosen parameters');
result=test_close(result,'lorentzfun imaginary',lor_i,[0 1/2],1e-15,1e-15,...
                  'zero-phase Lorentzian imaginary part is x/(1+x^2) for the chosen parameters');
result=test_close(result,'spden zero frequency',spden(2,10,0),1/300,1e-15,1e-15,...
                  'rank-two rotational spectral density at zero frequency is tau_c/(2L+1)');
result=test_close(result,'spden unit scaled frequency',spden(2,10,60),1/600,1e-15,1e-15,...
                  'when omega*tau_c is one, the Lorentzian spectral density is halved');

% Check minimum integer type selection at promotion boundaries
result=test_true(result,'min_int_type signed int8',strcmp(min_int_type(127,'signed'),'int8'),...
                 '127 is representable in signed 8-bit storage');
result=test_true(result,'min_int_type signed int16',strcmp(min_int_type(128,'signed'),'int16'),...
                 '128 requires signed 16-bit storage');
result=test_true(result,'min_int_type unsigned uint8',strcmp(min_int_type(255,'unsigned'),'uint8'),...
                 '255 is representable in unsigned 8-bit storage');
result=test_true(result,'min_int_type unsigned uint16',strcmp(min_int_type(256,'unsigned'),'uint16'),...
                 '256 requires unsigned 16-bit storage');

end


