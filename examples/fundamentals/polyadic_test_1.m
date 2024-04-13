% Unit tests for the polyadic object.
%
% i.kuprov@soton.ac.uk

function polyadic_test_1()

% Get random complex matrices
a=randn(7,7)+1i*randn(7,7);
b=randn(1,1)+1i*randn(1,1);
c=sprandn(9,9,0.1)+1i*sprandn(9,9,0.1);

% Get a random complex vector
v=randn(7*9,1)+1i*randn(7*9,1);

% Create the polyadic
P=polyadic({{polyadic({{a,b}}),c},{a b c}});

% Add prefixes and suffixes
p1=randn(size(P))+1i*randn(size(P));
p2=randn(size(P))+1i*randn(size(P));
s1=randn(size(P))+1i*randn(size(P));
s2=randn(size(P))+1i*randn(size(P));
P=prefix(p2,P); P=prefix(p1,P); 
P=suffix(P,s1); P=suffix(P,s2);

% Reference matrix
M=2*p1*p2*kron(a,kron(b,c))*s1*s2;

% Create-inflate test
if norm(M-inflate(P),1)<1e-3
    disp('Create-inflate test PASSED.');
else
    error('Create-inflate test FAILED.');
end

% Matrix-vector test
Pv=P*v; vP=v'*P; Mv=M*v; vM=v'*M;
N=norm(Pv-Mv,1)+norm(vP-vM,1);
if N<1e-3
    disp('Matrix-vector test PASSED.');
else
    error('Matrix-vector test FAILED.');
end

% Sparse test
if (issparse(inflate(P)))||(~issparse(full(P)))
    disp('Sparsity test PASSED.');
else
    error('Sparsity test FAILED.');
end

% Addition test 1
if norm(inflate(P-P),1)<1e-14
    disp('Addition test 1 PASSED.');
else
    error('Addition test 1 FAILED.');
end

% Addition test 2
if norm(inflate(2*P+P*3)-5*M,1)<1e-3
    disp('Addition test 2 PASSED.');
else
    error('Addition test 2 FAILED.');
end

% Conjugate-transpose test 1
if norm(P*v-(v'*P')')<1e-12
    disp('Conjugate-transpose test 1 PASSED.');
else
    error('Conjugate-transpose test 1 FAILED.');
end

% Conjugate-transpose test 2
if norm((inflate(P)')*v-inflate(P')*v)<1e-6
    disp('Conjugate-transpose test 2 PASSED.');
else
    error('Conjugate-transpose test 2 FAILED.');
end

% Size test
if norm(size(inflate(P))-size(P))<1e-12
    disp('Size test PASSED.');
else
    error('Size test FAILED.');
end

% Kron test 1
K=randn(2,2)+1i*randn(2,2);
if (norm(kron(inflate(P),K)-inflate(kron(P,K)),1)<1e-12)&&...
   (norm(kron(K,inflate(P))-inflate(kron(K,P)),1)<1e-12)
    disp('Kron test 1 PASSED.');
else
    error('Kron test 1 FAILED.');
end

% Kron test 2
K=polyadic({{randn(2,2)+1i*randn(2,2),randn(2,2)+1i*randn(2,2)}});
if (norm(kron(full(K),inflate(P))-inflate(kron(K,P)),1)<1e-12)&&...
   (norm(kron(inflate(P),full(K))-inflate(kron(P,K)),1)<1e-12)
    disp('Kron test 2 PASSED.');
else
    error('Kron test 2 FAILED.');
end

% Step test
spin_system=bootstrap('hush');
dt=1/norm(M,1);
v1=step(spin_system,P,v,dt);
v2=step(spin_system,M,v,dt);
if norm(v1-v2)<1e-12
    disp('Step test PASSED.');
else
    error('Step test FAILED.');
end

end

