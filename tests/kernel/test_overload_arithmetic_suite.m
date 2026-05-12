% Tests cheap overload arithmetic for cell, struct, RCV, and polyadic classes. Syntax:
%
%                    result=test_overload_arithmetic_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks elementwise cell overloads, recursive struct arithmetic,
% RCV sparse storage operations, and small polyadic arithmetic against
% explicit Matlab matrix references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_overload_arithmetic_suite()

% State the utility target of the test
result=new_test_result('kernel/overload_arithmetic_suite',...
                       'Cheap overload arithmetic',...
                       'cell, struct, RCV, and polyadic overloads must match explicit Matlab arithmetic on small examples.');

% Define small cell-array operands
A=[1 2;3 4];
B=[0 5;-1 2];
cell_a={A,B};
cell_b={eye(2),ones(2)};

% Check cell addition and subtraction overloads
cell_c=cell_a+cell_b;
result=test_close(result,'cell plus first element',cell_c{1},A+eye(2),1e-15,1e-15,...
                  'cell plus adds corresponding cells element by element');
result=test_close(result,'cell plus second element',cell_c{2},B+ones(2),1e-15,1e-15,...
                  'cell plus adds corresponding cells element by element');
cell_c=cell_a-cell_b;
result=test_close(result,'cell minus first element',cell_c{1},A-eye(2),1e-15,1e-15,...
                  'cell minus subtracts corresponding cells element by element');
cell_c=cell_a+2;
result=test_close(result,'cell scalar right plus',cell_c{1},cell_a{1}+2,1e-15,1e-15,...
                  'adding a scalar on the right applies it to every cell');
cell_c=2+cell_a;
result=test_close(result,'cell scalar left plus',cell_c{2},2+cell_a{2},1e-15,1e-15,...
                  'adding a scalar on the left applies it to every cell');

% Check cell scalar-array and matrix multiplication overloads
cell_c=cell_a.*[2 3];
result=test_close(result,'cell times right first element',cell_c{1},2*A,1e-15,1e-15,...
                  'cell times multiplies each cell by the matching scalar');
result=test_close(result,'cell times right second element',cell_c{2},3*B,1e-15,1e-15,...
                  'cell times multiplies each cell by the matching scalar');
cell_c=[2 3].*cell_a;
result=test_close(result,'cell times left second element',cell_c{2},3*B,1e-15,1e-15,...
                  'left scalar-array times multiplies each cell by the matching scalar');
R=diag([2 3]);
cell_c=cell_a*R;
result=test_close(result,'cell mtimes right',cell_c{1},A*R,1e-15,1e-15,...
                  'cell mtimes on the right multiplies every cell from the right');
cell_c=R*cell_a;
result=test_close(result,'cell mtimes left',cell_c{2},R*B,1e-15,1e-15,...
                  'cell mtimes on the left multiplies every cell from the left');

% Check cell totals and inflation shorthand
cell_sparse={sparse([1 0;0 2]),sparse([0 3;4 0])};
result=test_close(result,'cell totsum sparse',totsum(cell_sparse),sparse([1 3;4 2]),1e-15,1e-15,...
                  'totsum adds all cell entries and preserves sparse matrix arithmetic');
cell_inflated=inflate(cell_a);
result=test_close(result,'cell inflate first element',cell_inflated{1},A,1e-15,1e-15,...
                  'inflating a cell array of doubles leaves each numeric entry unchanged');
cell_complex=complex(cell_a);
result=test_close(result,'cell complex second element',cell_complex{2},complex(B),1e-15,1e-15,...
                  'complex applied to a cell array applies complex to every cell');

% Check recursive struct addition and left multiplication
s1.alpha=[1;2];
s1.beta.gamma=[3 4];
s2.alpha=[5;6];
s2.beta.gamma=[7 8];
s3=s1+s2;
result=test_close(result,'struct plus top field',s3.alpha,[6;8],1e-15,1e-15,...
                  'struct plus recursively adds matching numeric fields');
result=test_close(result,'struct plus nested field',s3.beta.gamma,[10 12],1e-15,1e-15,...
                  'struct plus recursively adds matching nested numeric fields');
s4=2*s1;
result=test_close(result,'struct mtimes top field',s4.alpha,[2;4],1e-15,1e-15,...
                  'numeric left multiplication is applied recursively to struct fields');
result=test_close(result,'struct mtimes nested field',s4.beta.gamma,[6 8],1e-15,1e-15,...
                  'numeric left multiplication is applied recursively to nested struct fields');

% Check RCV sparse storage and arithmetic against Matlab sparse references
S=sparse([1 0 2;0 3 0]);
T=sparse([0 4 0;5 0 6]);
R1=rcv(S);
R2=rcv(T);
result=test_close(result,'rcv full',full(R1),full(S),1e-15,1e-15,...
                  'full(rcv(S)) reconstructs the original sparse matrix');
result=test_close(result,'rcv plus rcv',full(R1+R2),full(S+T),1e-15,1e-15,...
                  'RCV plus concatenates entries and sparse conversion sums duplicates');
result=test_close(result,'rcv plus sparse',full(R1+T),full(S+T),1e-15,1e-15,...
                  'RCV plus Matlab sparse matches sparse addition');
result=test_close(result,'rcv minus',full(R1-R2),full(S-T),1e-15,1e-15,...
                  'RCV minus matches sparse subtraction');
result=test_close(result,'rcv scalar mtimes',full(2*R1),full(2*S),1e-15,1e-15,...
                  'scalar-matrix RCV multiplication scales stored values');
result=test_close(result,'rcv scalar times',full(3.*R1),full(3*S),1e-15,1e-15,...
                  'scalar-array RCV multiplication scales stored values');
result=test_close(result,'rcv rdivide',full(R1./2),full(S./2),1e-15,1e-15,...
                  'RCV right division by a scalar scales stored values');
result=test_close(result,'rcv transpose',full(R1.'),full(S.'),1e-15,1e-15,...
                  'RCV transpose swaps stored row and column indices');
C=sparse([1+1i 0;0 2-3i]);
result=test_close(result,'rcv ctranspose',full(rcv(C)'),full(C'),1e-15,1e-15,...
                  'RCV conjugate transpose swaps indices and conjugates values');
result=test_close(result,'rcv mtimes',full(R1*rcv(S.')),full(S*S.'),1e-15,1e-15,...
                  'RCV matrix multiplication matches Matlab sparse multiplication');
result=test_close(result,'rcv horzcat',full([R1 R2]),full([S T]),1e-15,1e-15,...
                  'horizontal RCV concatenation shifts right-hand column indices');
result=test_close(result,'rcv vertcat',full([R1;R2]),full([S;T]),1e-15,1e-15,...
                  'vertical RCV concatenation shifts lower row indices');
result=test_close(result,'rcv size vector',double(size(R1)),[2 3],1e-15,1e-15,...
                  'RCV size reports the stored matrix dimensions');

% Check small polyadic arithmetic against opened Kronecker products
P1=[1 2;0 3];
P2=[2 -1;4 1];
P3=[0 1;5 2];
P4=[3 0;-2 1];
P=polyadic({{P1,P2},{P3,P4}});
P_ref=kron(P1,P2)+kron(P3,P4);
result=test_close(result,'polyadic full',full(P),P_ref,1e-14,1e-14,...
                  'full(polyadic) opens and sums all stored Kronecker products');
result=test_close(result,'polyadic size',size(P),[4 4],1e-15,1e-15,...
                  'polyadic size is the product of core dimensions');
result=test_close(result,'polyadic nnz',nnz(P),nnz(P1)+nnz(P2)+nnz(P3)+nnz(P4),1e-15,1e-15,...
                  'polyadic nnz counts non-zero entries stored in all cores');
result=test_close(result,'polyadic plus matrix',full(P+speye(4)),P_ref+eye(4),1e-14,1e-14,...
                  'polyadic plus a matrix buffers a new additive term');
result=test_close(result,'polyadic scalar left mtimes',full(2*P),2*P_ref,1e-14,1e-14,...
                  'left scalar multiplication scales the buffered polyadic cores');
result=test_close(result,'polyadic scalar right mtimes',full(P*3),3*P_ref,1e-14,1e-14,...
                  'right scalar multiplication scales the buffered polyadic cores');
v=(1:4)';
result=test_close(result,'polyadic vector mtimes',P*v,P_ref*v,1e-14,1e-14,...
                  'polyadic multiplication by a dense vector matches the opened matrix');
result=test_close(result,'polyadic transpose',full(P.'),P_ref.',1e-14,1e-14,...
                  'polyadic transpose transposes every core and swaps prefixes with suffixes');
result=test_close(result,'polyadic ctranspose',full(P'),P_ref',1e-14,1e-14,...
                  'polyadic conjugate transpose conjugates and transposes every core');
result=test_close(result,'polyadic kron matrix',full(kron(polyadic({{P1}}),P2)),kron(P1,P2),1e-14,1e-14,...
                  'polyadic kron appends matrix cores without opening the product');

end


