function tests = test_prune_cpu
%TEST_PRUNE_CPU Unit tests for prune_cpu MEX.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Ensure compile
util_dir = fileparts(mfilename('fullpath'));
addpath(util_dir);

mexname = ['prune_cpu.' mexext];
if ~exist(fullfile(util_dir,mexname),'file')
    compile_prune_cpu();
end
end

function testQuantiseAndPruneComplexSparse(testCase)
% Large-ish sparse complex
rng(1);
m = 2000; n = 1500;
density = 2e-4;
A = sprand(m,n,density) + 1i*sprand(m,n,density);

% Add some values that should go to exact zeros after quantisation
A(10,10) = 0.49e-3 + 1i*0.49e-3;
A(20,20) = 0.49e-3 + 1i*0;

tol = 1e-3;

% Reference: quantise only stored entries, then explicitly drop zeros
[ii,jj,vv] = find(A);
vvq = tol * round((1/tol) * vv);
keep = (real(vvq)~=0) | (imag(vvq)~=0);
Aref = sparse(ii(keep), jj(keep), vvq(keep), m, n);

% Test target
Aout = prune_cpu(A, tol);

verifyTrue(testCase, issparse(Aout));
verifyTrue(testCase, ~isreal(Aout));

% Numerical equality (sparse equality ignores stored zeros)
verifyEqual(testCase, Aout, Aref);

% Pruning: nnz should match reference
verifyEqual(testCase, nnz(Aout), nnz(Aref));

% Explicitly check that the two injected near-zero entries disappeared
verifyEqual(testCase, full(Aout(10,10)), 0);
verifyEqual(testCase, full(Aout(20,20)), 0);
end

function testRealSparse(testCase)
rng(2);
A = sprand(1000,800,1e-3);
A(1,1) = 0.49e-2;
A(2,2) = 0.51e-2;

tol = 1e-2;
[ii,jj,vv] = find(A);
vvq = tol * round((1/tol) * vv);
keep = (vvq ~= 0);
Aref = sparse(ii(keep), jj(keep), vvq(keep), size(A,1), size(A,2));

Aout = prune_cpu(A, tol);

verifyTrue(testCase, issparse(Aout));
verifyTrue(testCase, isreal(Aout));
verifyEqual(testCase, Aout, Aref);
verifyEqual(testCase, nnz(Aout), nnz(Aref));
verifyEqual(testCase, full(Aout(1,1)), 0);
verifyEqual(testCase, full(Aout(2,2)), 0.01);
end

function testBadInputs(testCase)
Afull = rand(3);
tol = 1e-3;
verifyError(testCase, @() call_prune(Afull,tol), 'Spinach:prune_cpu:notSparse');
verifyError(testCase, @() call_prune(sparse(eye(3)), -1), 'Spinach:prune_cpu:tolVal');
verifyError(testCase, @() call_prune(sparse(eye(3)), [1 2]), 'Spinach:prune_cpu:tolType');

A = sparse(eye(3));
verifyError(testCase, @() call_prune_no_output(A,tol), 'Spinach:prune_cpu:nlhs');
end

function Aout = call_prune(A,tol)
Aout = prune_cpu(A,tol);
end

function call_prune_no_output(A,tol)
prune_cpu(A,tol);
end
