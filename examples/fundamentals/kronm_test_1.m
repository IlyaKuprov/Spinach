% Tests for the kron-times-matrix infrastructure.
%
% ilya.kuprov@weizmann.ac.il

function kronm_test_1()

% Pick random dimensions
nmats=randi(4,1,1)+2; disp(['number of matrices in Q: ' num2str(nmats)]);
ncols=randi(20,1,1); disp(['number of columns in x:  ' num2str(ncols)]);
dims=randi(7,nmats,1)+2; disp(['dims of matrices in Q:   ' num2str(dims')]);

% Generate real Q and x
Q_terms=cell(1,nmats); Q=1;
for n=1:nmats
    Q_terms{n}=randn(dims(n));
    Q=kron(Q,Q_terms{n});
end
x=randn(prod(dims),ncols);

% Compare the results
tic; x_a=kronm(Q_terms,x); disp(['kronm time:  ' num2str(toc()) ' seconds']);
tic; x_b=Q*x; disp(['mtimes time: ' num2str(toc()) ' seconds']);
if norm(x_a-x_b,1)<1e-6
    disp('Real matrix test PASSED.');
else
    error(['Real matrix test FAILED, norm(xa-xb)=' num2str(norm(x_a-x_b,1))]);
end

% Generate complex Q and x
Q_terms=cell(1,nmats); Q=1;
for n=1:nmats
    Q_terms{n}=randn(dims(n))+1i*randn(dims(n));
    Q=kron(Q,Q_terms{n});
end
x=randn(prod(dims),ncols)+1i*randn(prod(dims),ncols);

% Compare the results
tic; x_a=kronm(Q_terms,x); disp(['kronm time:  ' num2str(toc()) ' seconds']);
tic; x_b=Q*x; disp(['mtimes time: ' num2str(toc()) ' seconds']);
if norm(x_a-x_b,1)<1e-6
    disp('Complex matrix test PASSED.');
else
    error(['Complex matrix test FAILED, norm(xa-xb)=' num2str(norm(x_a-x_b,1))]);
end

end

