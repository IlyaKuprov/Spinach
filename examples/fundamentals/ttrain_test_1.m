% A simple test of ttclass object arithmetic.
%
% ilya.kuprov@weizmann.ac.il
% d.savostyanov@soton.ac.uk

function ttrain_test_1()

% Create a bunch of matrices
A=magic(5); A=A/norm(A,2);
B=randn(20); B=B/norm(B,2);
C=1i*rand(15); C=C/norm(C,2);

% Make their Kronecker product
P=ttclass(1,{A; B; C},0);

% Compute a function in TT
P_TT=full(P*P+3*P);

% Compute the usual way
P=kron(A,kron(B,C)); 
P_US=P*P+3*P;

% Compare results
if norm(P_TT-P_US,1)<100*eps('double')
    disp('Test passed.');
else
    error('Test failed.');
end

end

