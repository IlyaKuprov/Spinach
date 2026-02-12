% Delicate action and commutation tests for Hilbert-Liouville
% conversion and Stevens operators.
%
% ilya.kuprov@weizmann.ac.il

function commutation_8()

% Accuracy threshold
tol=1e-10;

% Test Hilbert-Liouville conversion identities
dim=6;
H=randn(dim)+1i*randn(dim);
R=randn(dim)+1i*randn(dim);
L_left=hilb2liouv(H,'left');
L_right=hilb2liouv(H,'right');
L_comm=hilb2liouv(H,'comm');
L_acomm=hilb2liouv(H,'acomm');
R_vec=hilb2liouv(R,'statevec');

% Build Hilbert-space products for vectorization checks
HR=H*R; RH=R*H;

% Compare left, right, commutator and anticommutator actions
err_left=norm(L_left*R_vec-HR(:),2);
err_right=norm(L_right*R_vec-RH(:),2);
err_comm=norm(L_comm*R_vec-reshape(HR-RH,[],1),2);
err_acomm=norm(L_acomm*R_vec-reshape(HR+RH,[],1),2);

% Report Hilbert-Liouville conversion failures
if max([err_left err_right err_comm err_acomm])>tol
    error('Hilbert-Liouville conversion test FAILED.');
end
disp('Hilbert-Liouville conversion test PASSED.');

% Test first-rank Stevens operators against Pauli operators
for mult=[2 3 5 8]

    % Build Stevens and Pauli operators
    L=pauli(mult);
    O_10=stevens(mult,1,0);
    O_11=stevens(mult,1,1);
    O_1m1=stevens(mult,1,-1);

    % Compare first-rank operators
    err_z=norm(O_10-L.z,'fro');
    err_x=norm(O_11-L.x,'fro');
    err_y=norm(O_1m1-L.y,'fro');

    % Report first-rank Stevens failures
    if max([err_z err_x err_y])>tol
        error('Stevens first-rank mapping test FAILED.');
    end

end
disp('Stevens first-rank mapping test PASSED.');

end

