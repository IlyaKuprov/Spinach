% Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax:
%
%          result=test_dynamic_overload_sparse_polyadic_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises object-operation dispatch for rcv, polyadic,
% and opium objects using small deterministic dense references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_overload_sparse_polyadic_suite()

% State the overload target of the test
result=new_test_result('kernel/dynamic_overload_sparse_polyadic_suite',...
                       'Dynamic sparse and polyadic overload dispatch',...
                       'RCV, polyadic, and OPIUM overloads must match explicit dense references on small deterministic objects.');

% Build deterministic sparse operands
A=sparse([1 0 2;0 -3 0]);
B=sparse([0 4 -1;5 0 6]);
ra=rcv(A);
rb=rcv(B);

% Exercise RCV construction, conversion, and size dispatch
result=test_close(result,'rcv sparse conversion',sparse(ra),A,1e-15,1e-15,...
                  'rcv sparse conversion must reproduce the input matrix');
result=test_close(result,'rcv full conversion',full(ra),full(A),1e-15,1e-15,...
                  'rcv full conversion must reproduce the dense input matrix');
result=test_close(result,'rcv size vector',double(size(ra)),double(size(A)),0,0,...
                  'rcv size must match the source sparse matrix');
result=test_close(result,'rcv first dimension',double(size(ra,1)),size(A,1),0,0,...
                  'rcv size along the first dimension must match Matlab');

% Exercise RCV arithmetic dispatch
result=test_close(result,'rcv plus rcv',full(ra+rb),full(A+B),1e-15,1e-15,...
                  'rcv plus must add two RCV sparse matrices');
result=test_close(result,'rcv plus sparse',full(ra+B),full(A+B),1e-15,1e-15,...
                  'rcv plus must add an RCV sparse matrix and a Matlab sparse matrix');
result=test_close(result,'rcv minus',full(ra-rb),full(A-B),1e-15,1e-15,...
                  'rcv minus must subtract two RCV sparse matrices');
result=test_close(result,'rcv array times',full(3.*ra),full(3*A),1e-15,1e-15,...
                  'rcv times must scale stored values by a scalar');
result=test_close(result,'rcv right mtimes scalar',full(ra*2),full(A*2),1e-15,1e-15,...
                  'rcv mtimes must scale by a right scalar');
result=test_close(result,'rcv left mtimes scalar',full(2*ra),full(2*A),1e-15,1e-15,...
                  'rcv mtimes must scale by a left scalar');
result=test_close(result,'rcv rdivide scalar',full(ra./2),full(A/2),1e-15,1e-15,...
                  'rcv rdivide must divide stored values by a scalar');

% Exercise direct RCV method-name dispatch for coverage tracking
result=test_close(result,'rcv direct plus',full(plus(ra,rb)),full(A+B),1e-15,1e-15,...
                  'direct plus() must dispatch to the RCV overload');
result=test_close(result,'rcv direct minus',full(minus(ra,rb)),full(A-B),1e-15,1e-15,...
                  'direct minus() must dispatch to the RCV overload');
result=test_close(result,'rcv direct times',full(times(3,ra)),full(3*A),1e-15,1e-15,...
                  'direct times() must dispatch to the RCV overload');
result=test_close(result,'rcv direct mtimes',full(mtimes(ra,2)),full(A*2),1e-15,1e-15,...
                  'direct mtimes() must dispatch to the RCV overload');
result=test_close(result,'rcv direct rdivide',full(rdivide(ra,2)),full(A/2),1e-15,1e-15,...
                  'direct rdivide() must dispatch to the RCV overload');

% Exercise RCV matrix product and transposition dispatch
C=sparse([2 0;-1 3;0 4]);
L=sparse([2 0;-1 3;0 4]);
result=test_close(result,'rcv mtimes sparse right',ra*C,full(A*C),1e-15,1e-15,...
                  'rcv mtimes must multiply a Matlab sparse matrix on the right');
result=test_close(result,'rcv mtimes sparse left',L*ra,full(L*A),1e-15,1e-15,...
                  'rcv mtimes must multiply a Matlab sparse matrix on the left');
result=test_close(result,'rcv transpose',full(ra.'),full(A.'),1e-15,1e-15,...
                  'rcv transpose must swap rows and columns');
complex_rcv=rcv(A+1i*B);
result=test_close(result,'rcv ctranspose',full(complex_rcv'),full((A+1i*B)'),1e-15,1e-15,...
                  'rcv ctranspose must conjugate values and swap rows and columns');
result=test_close(result,'rcv horizontal concatenation',full([ra rb]),full([A B]),1e-15,1e-15,...
                  'rcv horzcat must concatenate columns');
result=test_close(result,'rcv vertical concatenation',full([ra;rb]),full([A;B]),1e-15,1e-15,...
                  'rcv vertcat must concatenate rows');
result=test_close(result,'rcv direct ctranspose',full(ctranspose(complex_rcv)),...
                  full((A+1i*B)'),1e-15,1e-15,...
                  'direct ctranspose() must dispatch to the RCV overload');
result=test_close(result,'rcv direct horzcat',full(horzcat(ra,rb)),full([A B]),1e-15,1e-15,...
                  'direct horzcat() must dispatch to the RCV overload');
result=test_close(result,'rcv direct vertcat',full(vertcat(ra,rb)),full([A;B]),1e-15,1e-15,...
                  'direct vertcat() must dispatch to the RCV overload');
result=test_close(result,'rcv gather cpu no-op',full(gather(ra)),full(A),1e-15,1e-15,...
                  'gather on a CPU RCV object must be a no-op');

% Exercise RCV sparsity plotting dispatch offscreen
set(0,'DefaultFigureVisible','off');
figures_before=numel(findall(0,'Type','figure'));
spy(ra);
figures_after=numel(findall(0,'Type','figure'));
close all;
result=test_true(result,'rcv direct spy',figures_after>=figures_before+1,...
                 'direct spy() must dispatch to the RCV overload and create a figure');

% Exercise guarded RCV GPU dispatch when a GPU is available
if local_have_gpu()
    gpu_rcv=gpuArray(ra);
    if ~gpu_rcv.isGPU
        error('FAILED: rcv gpuArray did not mark the object as GPU-resident.');
    end
    cpu_rcv=gather(gpu_rcv);
    result=test_close(result,'rcv gpuArray gather round-trip',full(cpu_rcv),full(A),1e-15,1e-15,...
                      'rcv GPU upload followed by gather must preserve values');
else
    result.messages{end+1}='SKIP: rcv gpuArray path skipped because no usable GPU is available.';
end

% Build deterministic polyadic operands
P=polyadic({{[1 2;0 -1],[0 3;4 5]},...
            {[2 0;1 1],[-1 2;3 0]}});
Q=polyadic({{[0 1;-2 3],[5 -1;0 2]}});
p_ref=kron([1 2;0 -1],[0 3;4 5])+...
      kron([2 0;1 1],[-1 2;3 0]);
q_ref=kron([0 1;-2 3],[5 -1;0 2]);

% Exercise polyadic construction, conversion, and structural predicates
validate(P);
result.messages{end+1}='PASS: polyadic validate accepted the deterministic object.';
result=test_close(result,'polyadic full',full(P),p_ref,1e-15,1e-15,...
                  'polyadic full must open the buffered Kronecker sum');
result=test_close(result,'polyadic size vector',double(size(P)),double(size(p_ref)),0,0,...
                  'polyadic size must match the opened dense matrix');
if isempty(P)||~allfinite(P)||nnz(P)~=12
    error('FAILED: polyadic structural predicates returned unexpected values.');
end
result.messages{end+1}='PASS: polyadic isempty, allfinite, and nnz predicates returned expected values.';

% Exercise polyadic arithmetic and Kronecker dispatch
result=test_close(result,'polyadic plus',full(P+Q),p_ref+q_ref,1e-15,1e-15,...
                  'polyadic plus must buffer sums without changing values');
result=test_close(result,'polyadic minus',full(P-Q),p_ref-q_ref,1e-15,1e-15,...
                  'polyadic minus must buffer signed sums without changing values');
result=test_close(result,'polyadic scalar left mtimes',full(2*P),2*p_ref,1e-15,1e-15,...
                  'left scalar multiplication must scale a polyadic');
result=test_close(result,'polyadic scalar right mtimes',full(P*3),3*p_ref,1e-15,1e-15,...
                  'right scalar multiplication must scale a polyadic');
vec=[1;-1;2;0];
result=test_close(result,'polyadic dense-vector mtimes',P*vec,p_ref*vec,1e-15,1e-15,...
                  'polyadic mtimes must act on dense vectors');
result=test_close(result,'polyadic direct plus',full(plus(P,Q)),p_ref+q_ref,1e-15,1e-15,...
                  'direct plus() must dispatch to the polyadic overload');
result=test_close(result,'polyadic direct minus',full(minus(P,Q)),p_ref-q_ref,1e-15,1e-15,...
                  'direct minus() must dispatch to the polyadic overload');
result=test_close(result,'polyadic direct mtimes',mtimes(P,vec),p_ref*vec,1e-15,1e-15,...
                  'direct mtimes() must dispatch to the polyadic overload');
result=test_close(result,'polyadic sparse prefix mtimes',full(sparse(eye(4))*P),p_ref,1e-15,1e-15,...
                  'sparse left multiplication must attach a value-preserving prefix');
result=test_close(result,'polyadic sparse suffix mtimes',full(P*sparse(eye(4))),p_ref,1e-15,1e-15,...
                  'sparse right multiplication must attach a value-preserving suffix');
result=test_close(result,'polyadic kron right numeric',full(kron(P,[1 2;3 4])),...
                  kron(p_ref,[1 2;3 4]),1e-15,1e-15,...
                  'polyadic kron with a numeric right operand must match dense kron');
result=test_close(result,'polyadic kron left numeric',full(kron([1 0;0 -1],Q)),...
                  kron([1 0;0 -1],q_ref),1e-15,1e-15,...
                  'polyadic kron with a numeric left operand must match dense kron');

% Exercise polyadic transpose, prefix, suffix, simplify, and inflate dispatch
result=test_close(result,'polyadic transpose',full(P.'),p_ref.',1e-15,1e-15,...
                  'polyadic transpose must transpose every represented matrix factor');
complex_poly=polyadic({{[1 1i;0 2],[3 0;-1i 4]}});
result=test_close(result,'polyadic ctranspose',full(complex_poly'),full(complex_poly)',1e-15,1e-15,...
                  'polyadic ctranspose must conjugate-transpose the represented matrix');
result=test_close(result,'polyadic direct ctranspose',full(ctranspose(complex_poly)),...
                  full(complex_poly)',1e-15,1e-15,...
                  'direct ctranspose() must dispatch to the polyadic overload');
L=sparse([1 0 0 0;0 2 0 0;0 0 3 0;0 0 0 4]);
R=sparse([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]);
result=test_close(result,'polyadic prefix',full(prefix(L,P)),full(L*p_ref),1e-15,1e-15,...
                  'polyadic prefix must multiply represented values on the left');
result=test_close(result,'polyadic suffix',full(suffix(P,R)),full(p_ref*R),1e-15,1e-15,...
                  'polyadic suffix must multiply represented values on the right');
result=test_close(result,'polyadic simplify',full(simplify(prefix(speye(4),suffix(P,speye(4))))),...
                  p_ref,1e-15,1e-15,...
                  'polyadic simplify must remove identity prefixes and suffixes');
result=test_close(result,'polyadic inflate',inflate(P),p_ref,1e-15,1e-15,...
                  'polyadic inflate must produce the opened dense representation');

% Exercise guarded polyadic GPU dispatch when a GPU is available
if local_have_gpu()
    gpu_poly=gpuArray(P);
    result=test_close(result,'polyadic gpuArray full gather',gather(full(gpu_poly)),...
                      p_ref,1e-15,1e-15,...
                      'polyadic GPU upload followed by full and gather must preserve values');
else
    result.messages{end+1}='SKIP: polyadic gpuArray path skipped because no usable GPU is available.';
end

% Exercise OPIUM construction, conversion, size, multiplication, and kron
unit_op=opium(3,1);
scaled_op=opium(3,2);
other_op=opium(3,-4);
M=[1 2 3;4 5 6;7 8 9];
result=test_close(result,'opium size vector',double(size(scaled_op)),[3 3],0,0,...
                  'opium size must report a square matrix of its stored dimension');
result=test_close(result,'opium full unit',full(unit_op),eye(3),1e-15,1e-15,...
                  'opium full must open a unit OPIUM as an identity matrix');
result=test_close(result,'opium sparse scaled',sparse(scaled_op),2*speye(3),1e-15,1e-15,...
                  'opium sparse must include the stored coefficient');
result=test_close(result,'opium left scalar mtimes',sparse(5*scaled_op),10*speye(3),1e-15,1e-15,...
                  'left scalar multiplication must scale the OPIUM coefficient');
result=test_close(result,'opium right scalar mtimes',sparse(scaled_op*5),10*speye(3),1e-15,1e-15,...
                  'right scalar multiplication must scale the OPIUM coefficient');
result=test_close(result,'opium dense left mtimes',M*scaled_op,2*M,1e-15,1e-15,...
                  'dense left multiplication by OPIUM must scale the dense matrix');
result=test_close(result,'opium dense right mtimes',scaled_op*M,2*M,1e-15,1e-15,...
                  'dense right multiplication by OPIUM must scale the dense matrix');
result=test_close(result,'opium opium mtimes',sparse(scaled_op*other_op),-8*speye(3),1e-15,1e-15,...
                  'OPIUM-by-OPIUM multiplication must multiply coefficients');
result=test_close(result,'opium kron opium',sparse(kron(scaled_op,opium(2,3))),...
                  6*speye(6),1e-15,1e-15,...
                  'OPIUM kron OPIUM must multiply dimensions and coefficients');
result=test_close(result,'opium kron numeric right',kron(scaled_op,[1 2;3 4]),...
                  kron(2*speye(3),[1 2;3 4]),1e-15,1e-15,...
                  'OPIUM kron numeric must match dense Kronecker expansion');
if ~allfinite(scaled_op)||~isnumeric(scaled_op)||~ismatrix(scaled_op)||nnz(scaled_op)~=1
    error('FAILED: opium scalar predicates returned unexpected values.');
end
result.messages{end+1}='PASS: opium allfinite, isnumeric, ismatrix, and nnz predicates returned expected values.';

end


function answer=local_have_gpu()

% Default to no GPU
answer=false();

% Leave silently when the GPU API is absent
if exist('gpuDevice','file')~=2
    return
end

% Try selecting the current GPU
try
    gpuDevice();
    answer=true();
catch
    answer=false();
end

end

