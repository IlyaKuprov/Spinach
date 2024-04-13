% Internal consistency of Blicharski invariants and their 
% polarisation relationships with inner products of irre-
% ducible spherical tensor coefficients, second rank.
%
% i.kuprov@soton.ac.uk

function blinv_test()

% Random tensors
A=rand(3,3); B=rand(3,3);

% Spherical expansions
[~,~,Phi_A]=mat2sphten(A);
[~,~,Phi_B]=mat2sphten(B);

% Blicharski invariants
[~,Dsq_A]=blinv(A);
[~,Dsq_B]=blinv(B);
[~,X_AB]=blprod(A,B);

% Blicharski against Phi products, rank 2 self
DA=(Phi_A(3)^2-2*Phi_A(2)*Phi_A(4)...
              +2*Phi_A(1)*Phi_A(5))-(2/3)*Dsq_A;
DB=(Phi_B(3)^2-2*Phi_B(2)*Phi_B(4)...
              +2*Phi_B(1)*Phi_B(5))-(2/3)*Dsq_B;
if (abs(DA)>10*eps)||(abs(DB)>10*eps)
    error('Test 1 failed.');
else
    disp('Test 1 passed.');
end

% Blicharski against Phi product, rank 2 cross
D_AB=(Phi_A(1)*Phi_B(5)-Phi_A(2)*Phi_B(4)+...
      Phi_A(3)*Phi_B(3)-Phi_A(4)*Phi_B(2)+...
      Phi_A(5)*Phi_B(1))-(2/3)*X_AB;
if (abs(D_AB)>10*eps)
    error('Test 2 failed.');
else
    disp('Test 2 passed.');
end

end

