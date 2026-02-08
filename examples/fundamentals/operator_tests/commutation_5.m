% Tests of basic product and commutation relations between
% bosonic operators.
%
% ilya.kuprov@weizmann.ac.il

function commutation_5()

% Random dimension
nlevels=5+randi(20);

% Accuracy threshold
acc=10*nlevels*eps('double');

% Weyl operators
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
    error('Bosonic operator commutation test 1 FAILED.');
else
    disp('Bosonic operator commutation test 1 PASSED.');
end

% Boundary commutator test
C=comm(A.a,A.c); C(end,end)=1;
err=norm(C-A.u,'fro')/norm(A.u,'fro');
if err>acc
    error('Bosonic operator commutation test 2 FAILED.');
else
    disp('Bosonic operator commutation test 2 PASSED.');
end

% Bosonic monomials
B=boson_mono(nlevels);

% Exact commutator test
for n=1:numel(B)

    % Physical indices
    [k,q]=lin2kq(nlevels,n);

    % Deviation norm testing
    err=norm(comm(A.n,B{n})-(k-q)*B{n},'fro')/norm(B{n},'fro');
    if err>acc
        error('Bosonic operator commutation test 3 FAILED.');
    end

end
disp('Bosonic operator commutation test 3 PASSED.');

end

