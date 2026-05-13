% Tests dynamic dispatch of tensor-train class overloads. Syntax:
%
%              result=test_dynamic_overload_ttclass_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises ttclass object operations against exact dense
% references on one-core and two-core tensor trains.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_overload_ttclass_suite()

% Announce the test target
fprintf('TESTING: Tensor-train overload dispatch\n');

% State the overload target of the test
result=new_test_result('kernel/dynamic_overload_ttclass_suite',...
                       'Dynamic ttclass overload dispatch',...
                       'ttclass overloads must match explicit dense references on small deterministic tensor trains.');

% Build one-core tensor trains and dense references
A=[1 2;3 4];
B=[0 5;6 1];
T=ttclass(2,{A},0);
U=ttclass(-3,{B},0);
t_ref=2*A;
u_ref=-3*B;

% Exercise constructor, full, size, sizes, ranks, numel, and subsref dispatch
result=test_close(result,'ttclass full one-core',full(T),t_ref,1e-15,1e-15,...
                  'ttclass full must open a one-core tensor train');
result=test_close(result,'ttclass size vector',double(size(T)),double(size(t_ref)),0,0,...
                  'ttclass size must match the dense matrix size');
result=test_close(result,'ttclass first dimension',size(T,1),size(t_ref,1),0,0,...
                  'ttclass size along the first dimension must match Matlab');
result=test_close(result,'ttclass mode sizes',T.sizes,[2 2],0,0,...
                  'ttclass sizes must report one physical row and column mode');
result=test_close(result,'ttclass ranks',T.ranks,[1;1],0,0,...
                  'ttclass ranks must report unit boundary ranks for a rank-one train');
result=test_close(result,'ttclass direct sizes',sizes(T),[2 2],0,0,...
                  'direct sizes() must dispatch to the ttclass overload');
result=test_close(result,'ttclass direct ranks',ranks(T),[1;1],0,0,...
                  'direct ranks() must dispatch to the ttclass overload');
result=test_close(result,'ttclass numel',numel(T),numel(t_ref),0,0,...
                  'ttclass numel must match the represented matrix element count');
result=test_close(result,'ttclass scalar subsref',T(2,1),t_ref(2,1),1e-15,1e-15,...
                  'ttclass scalar indexing must match dense indexing');
idx=substruct('()',{2,1});
result=test_close(result,'ttclass direct subsref',subsref(T,idx),t_ref(2,1),1e-15,1e-15,...
                  'direct subsref() must dispatch to the ttclass overload');
if (T.ncores~=1)||(T.ntrains~=1)||~ismatrix(T)||~isnumeric(T)||~isreal(T)
    error('FAILED: ttclass structural predicates returned unexpected values.');
end
result.messages{end+1}='PASS: ttclass ncores, ntrains, ismatrix, isnumeric, and isreal predicates returned expected values.';

% Exercise addition, subtraction, scalar multiplication, and scalar division
result=test_close(result,'ttclass plus',full(T+U),t_ref+u_ref,1e-15,1e-15,...
                  'ttclass plus must concatenate buffered trains without changing values');
result=test_close(result,'ttclass minus',full(T-U),t_ref-u_ref,1e-15,1e-15,...
                  'ttclass minus must concatenate signed buffered trains without changing values');
result=test_close(result,'ttclass scalar left mtimes',full(4*T),4*t_ref,1e-15,1e-15,...
                  'left scalar multiplication must scale a tensor train');
result=test_close(result,'ttclass scalar right mtimes',full(T*4),4*t_ref,1e-15,1e-15,...
                  'right scalar multiplication must scale a tensor train');
result=test_close(result,'ttclass rdivide scalar',full(T./2),t_ref/2,1e-15,1e-15,...
                  'element-wise scalar division must scale a tensor train');
result=test_close(result,'ttclass mrdivide scalar',full(T/2),t_ref/2,1e-15,1e-15,...
                  'matrix scalar division must scale a tensor train');
result=test_close(result,'ttclass direct plus',full(plus(T,U)),t_ref+u_ref,1e-15,1e-15,...
                  'direct plus() must dispatch to the ttclass overload');
result=test_close(result,'ttclass direct minus',full(minus(T,U)),t_ref-u_ref,1e-15,1e-15,...
                  'direct minus() must dispatch to the ttclass overload');
result=test_close(result,'ttclass direct rdivide',full(rdivide(T,2)),t_ref/2,1e-15,1e-15,...
                  'direct rdivide() must dispatch to the ttclass overload');
result=test_close(result,'ttclass direct mrdivide',full(mrdivide(T,2)),t_ref/2,1e-15,1e-15,...
                  'direct mrdivide() must dispatch to the ttclass overload');

% Exercise matrix, vector, and tensor-train multiplication dispatch
rhs_vec=[1;-2];
result=test_close(result,'ttclass dense-vector mtimes',T*rhs_vec,t_ref*rhs_vec,1e-15,1e-15,...
                  'ttclass mtimes must act on dense vectors');
result=test_close(result,'ttclass dense-matrix mtimes',T*[1 0;0 -1],t_ref*[1 0;0 -1],1e-15,1e-15,...
                  'ttclass mtimes must act on dense matrices');
result=test_close(result,'ttclass tensor mtimes',full(T*U),t_ref*u_ref,1e-12,1e-12,...
                  'ttclass mtimes must multiply two tensor trains');
result=test_close(result,'ttclass direct mtimes',full(mtimes(T,U)),t_ref*u_ref,1e-12,1e-12,...
                  'direct mtimes() must dispatch to the ttclass overload');
result=test_close(result,'ttclass dot product object',full(dot(T,U)),t_ref'*u_ref,1e-12,1e-12,...
                  'ttclass dot must dispatch through Hermitian tensor-train multiplication');
result=test_close(result,'ttclass Hadamard dot',hdot(T,U),sum(conj(t_ref(:)).*u_ref(:)),1e-12,1e-12,...
                  'ttclass hdot must match the dense Frobenius inner product');

% Exercise conjugation, transposition, trace, diagonal, sum, and mean dispatch
Z=ttclass(1+2i,{[1 2i;-3i 4]},0);
z_ref=full(Z);
result=test_close(result,'ttclass conj',full(conj(Z)),conj(z_ref),1e-15,1e-15,...
                  'ttclass conj must conjugate all cores and coefficients');
result=test_close(result,'ttclass transpose',full(Z.'),z_ref.',1e-15,1e-15,...
                  'ttclass transpose must swap physical matrix indices');
result=test_close(result,'ttclass ctranspose',full(Z'),z_ref',1e-15,1e-15,...
                  'ttclass ctranspose must conjugate and swap physical matrix indices');
result=test_close(result,'ttclass direct ctranspose',full(ctranspose(Z)),z_ref',1e-15,1e-15,...
                  'direct ctranspose() must dispatch to the ttclass overload');
if isreal(Z)
    error('FAILED: ttclass isreal returned true for a complex tensor train.');
end
result.messages{end+1}='PASS: ttclass isreal detects complex tensor trains.';
result=test_close(result,'ttclass trace',trace(T),trace(t_ref),1e-15,1e-15,...
                  'ttclass trace must match the dense trace');
result=test_close(result,'ttclass diag matrix to vector',full(diag(T)),diag(t_ref),1e-15,1e-15,...
                  'ttclass diag must extract the represented matrix diagonal');
tt_vec=ttclass(1,{[2;5]},0);
result=test_close(result,'ttclass diag vector to matrix',full(diag(tt_vec)),diag([2;5]),1e-15,1e-15,...
                  'ttclass diag must build a diagonal tensor-train matrix from a vector');
result=test_close(result,'ttclass sum dim one',full(sum(T,1)),sum(t_ref,1),1e-15,1e-15,...
                  'ttclass sum along rows must match dense summation');
result=test_close(result,'ttclass sum dim two',full(sum(T,2)),sum(t_ref,2),1e-15,1e-15,...
                  'ttclass sum along columns must match dense summation');
result=test_close(result,'ttclass mean dim one',full(mean(T,1)),mean(t_ref,1),1e-15,1e-15,...
                  'ttclass mean along rows must match dense means');
result=test_close(result,'ttclass mean dim two',full(mean(T,2)),mean(t_ref,2),1e-15,1e-15,...
                  'ttclass mean along columns must match dense means');

% Exercise coefficient clearing, packing, orthogonalisation, truncation, and shrinkage
cleared=clearcoeff(T);
result=test_close(result,'ttclass clearcoeff value',full(cleared),t_ref,1e-15,1e-15,...
                  'clearcoeff must preserve the represented matrix');
result=test_close(result,'ttclass clearcoeff coefficient',cleared.coeff,1,1e-15,1e-15,...
                  'clearcoeff must absorb the physical coefficient into the core');
packed=pack(T+U);
result=test_close(result,'ttclass pack',full(packed),t_ref+u_ref,1e-12,1e-12,...
                  'pack must absorb the addition buffer without changing values');
ortho=ttort(T,+1);
result=test_close(result,'ttclass ttort value',full(ortho),t_ref,1e-12,1e-12,...
                  'ttort must preserve the represented value');
[normalised,lognrm]=ttort(T,+1);
result=test_close(result,'ttclass ttort logarithmic norm',exp(lognrm)*full(normalised),...
                  t_ref,1e-12,1e-12,...
                  'ttort with a log-norm output must preserve values after rescaling');
if ~isfinite(lognrm)
    error('FAILED: ttclass ttort returned a non-finite log norm.');
end
result.messages{end+1}='PASS: ttclass ttort returned a finite logarithmic norm.';
truncated=truncate(ortho);
result=test_close(result,'ttclass truncate value',full(truncated),t_ref,1e-12,1e-12,...
                  'truncate must preserve an exactly representable one-core train');
shrunk=shrink(T+U);
result=test_close(result,'ttclass shrink value',full(shrunk),t_ref+u_ref,1e-12,1e-12,...
                  'shrink must compress a one-core buffered sum without changing values');
result=test_close(result,'ttclass Frobenius norm',norm(T,'fro'),norm(t_ref,'fro'),1e-12,1e-12,...
                  'ttclass Frobenius norm must match the dense Frobenius norm');

% Exercise direct AMEn summation and isolated SPMD-save wrapper dispatch
sum_opts=struct('max_swp',10,'init_guess_rank',1,'enrichment_rank',0,'verb',0);
sum_train=ttclass(2,{1;1},0)+ttclass(3,{1;1},0);
summed=amensum(sum_train,1e-12,sum_opts);
result=test_close(result,'ttclass direct amensum',full(summed),5,1e-12,1e-12,...
                  'direct amensum() must compress a buffered two-core scalar sum');
save_file=[tempname '.mat'];
save_anyway(save_file,t_ref);
saved_data=load(save_file,'variable');
delete(save_file);
result=test_close(result,'save_anyway round trip',saved_data.variable,t_ref,1e-15,1e-15,...
                  'save_anyway() must write the supplied variable under the documented name');

% Exercise unit_like, rand, kron, and vec on object instances
result=test_close(result,'ttclass unit_like',full(unit_like(T)),eye(2),1e-15,1e-15,...
                  'unit_like must build an identity tensor train with matching topology');
rng(20240512);
rand_train=rand(ttclass(1,{A;B},0),2);
if ~isa(rand_train,'ttclass')||any(rand_train.ranks~=[1;2;1])||~all(isfinite(full(rand_train)),'all')
    error('FAILED: ttclass rand did not produce a finite rank-two tensor train.');
end
result.messages{end+1}='PASS: ttclass rand produced a finite rank-two tensor train with the requested rank.';
result=test_close(result,'ttclass kron one-core',full(kron(T,U)),kron(t_ref,u_ref),1e-12,1e-12,...
                  'one-core ttclass kron must match dense Kronecker multiplication');
result=test_close(result,'ttclass vec one-core',full(vec(T)),t_ref(:),1e-15,1e-15,...
                  'one-core ttclass vec must match dense column vectorisation');

% Exercise two-core full, multiplication, vectorisation, and bit-reversal paths
C=[2 -1;0 3];
D=[1 0;0 -2];
T2=ttclass(2,{A;B},0);
U2=ttclass(-1,{C;D},0);
two_ref=2*kron(A,B);
u_two_ref=-kron(C,D);
result=test_close(result,'ttclass full two-core',full(T2),two_ref,1e-15,1e-15,...
                  'ttclass full must open a two-core tensor train');
result=test_close(result,'ttclass two-core mtimes',full(T2*U2),two_ref*u_two_ref,1e-10,1e-10,...
                  'ttclass multiplication must handle two-core tensor trains');
result=test_close(result,'ttclass vec two-core',full(vec(T2)),2*kron(A(:),B(:)),1e-15,1e-15,...
                  'two-core ttclass vec must follow tensor-train Kronecker vectorisation order');
result=test_close(result,'ttclass revert two-core',full(revert(T2)),2*kron(B,A),1e-15,1e-15,...
                  'ttclass revert must reverse the tensor-train core order');

% Exercise AMEn solve on the smallest deterministic tensor train system
lhs=ttclass(2,{1},0);
rhs=ttclass(6,{1},0);
init=ttclass(1,{1},0);
opts=struct('nswp',2,'verb',0,'enrichment_rank',0,'max_full_size',10);
solution=amensolve(lhs,rhs,1e-12,opts,init);
result=test_close(result,'ttclass amensolve scalar system',full(solution),3,1e-10,1e-10,...
                  'AMEn solve must recover the exact scalar solution for a one-core system');
try
    mldivide(T,3);
    error('FAILED: direct mldivide() accepted a non-tensor right-hand side.');
catch err
    result=test_true(result,'ttclass direct mldivide grumbler',...
                     contains(err.message,'both arguments should be tensor trains'),...
                     'direct mldivide() must reject non-tensor right-hand sides');
end

end
