% Tests of basic product and commutation relations between
% bosonic operators.
%
% ilya.kuprov@weizmann.ac.il

function commutation_5()

% Random dimension
nlevels=5+randi(20);

% Accuracy threshold
acc=10*nlevels*eps('double');

% Operators
A=weyl(nlevels);

% Product test
err=norm(A.c*A.a-A.n,'fro')/norm(A.n,'fro');
if err>acc
    error('Bosonic operator product test FAILED.');
else
    disp('Bosonic operator product test PASSED.');
end

% Exact commutator test
err_a=norm(comm(A.n,A.c)-A.c,'fro')/norm(A.c,'fro');
err_b=norm(comm(A.n,A.a)+A.a,'fro')/norm(A.a,'fro');
if (err_a>acc)||(err_b>acc)
    error('Bosonic operator commutation test FAILED.');
else
    disp('Bosonic operator commutation test PASSED.');
end

% Boundary commutator test
C=comm(A.a,A.c); C(end,end)=1;
err=norm(C-A.u,'fro')/norm(A.u,'fro');
if err>acc
    error('Bosonic operator commutation test FAILED.');
else
    disp('Bosonic operator commutation test PASSED.');
end

end

