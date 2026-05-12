% Tests dynamic dispatch of cheap cell, struct, and double overloads. Syntax:
%
%              result=test_dynamic_overload_cell_struct_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises object-operation dispatch for cell arithmetic,
% cell multiplication, cell utility overloads, recursive structure
% arithmetic, and the double inflate no-op.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_overload_cell_struct_suite()

% State the overload target of the test
result=new_test_result('tmp/dynamic_overloads_retry/cell_struct',...
                       'Dynamic cell and struct overload dispatch',...
                       'Cell, struct, and double overloads must match explicit dense references on small deterministic objects.');

% Build deterministic matrix operands
A=[1 2;3 4];
B=[0 5;-1 2];
C=[2 -3;4 1];
D=[-2 0;5 3];
cell_a={A,B};
cell_b={C,D};

% Exercise cell plus and minus through operator dispatch
cell_c=cell_a+cell_b;
result=test_close(result,'cell plus first element',cell_c{1},A+C,1e-15,1e-15,...
                  'cell plus must add corresponding cells');
result=test_close(result,'cell plus second element',cell_c{2},B+D,1e-15,1e-15,...
                  'cell plus must add corresponding cells');
cell_c=cell_a-cell_b;
result=test_close(result,'cell minus first element',cell_c{1},A-C,1e-15,1e-15,...
                  'cell minus must subtract corresponding cells');
result=test_close(result,'cell minus second element',cell_c{2},B-D,1e-15,1e-15,...
                  'cell minus must subtract corresponding cells');

% Exercise numeric-cell plus and minus dispatch
cell_c=10+cell_a;
result=test_close(result,'cell left numeric plus',cell_c{1},10+A,1e-15,1e-15,...
                  'numeric left plus must add to each cell');
cell_c=cell_a-3;
result=test_close(result,'cell right numeric minus',cell_c{2},B-3,1e-15,1e-15,...
                  'numeric right minus must subtract from each cell');
cell_c=7-cell_a;
result=test_close(result,'cell left numeric minus',cell_c{1},7-A,1e-15,1e-15,...
                  'numeric left minus must subtract each cell from the scalar');

% Exercise cell times and mtimes dispatch
weights=[2 3];
cell_c=cell_a.*weights;
result=test_close(result,'cell right array times',cell_c{2},3*B,1e-15,1e-15,...
                  'right array times must scale matching cells');
cell_c=weights.*cell_a;
result=test_close(result,'cell left array times',cell_c{1},2*A,1e-15,1e-15,...
                  'left array times must scale matching cells');
R=[2 1;0 -1];
cell_c=cell_a*R;
result=test_close(result,'cell right mtimes',cell_c{1},A*R,1e-15,1e-15,...
                  'cell mtimes must right-multiply every cell');
cell_c=R*cell_a;
result=test_close(result,'cell left mtimes',cell_c{2},R*B,1e-15,1e-15,...
                  'cell mtimes must left-multiply every cell');

% Exercise cell utility overloads through real objects
cell_sparse={sparse(A),sparse(B)};
result=test_close(result,'cell totsum',totsum(cell_sparse),sparse(A+B),1e-15,1e-15,...
                  'totsum must add all cell entries');
cell_complex=complex(cell_a);
result=test_close(result,'cell complex',cell_complex{1},complex(A),1e-15,1e-15,...
                  'complex must apply to each cell entry');
cell_inflated=inflate(cell_a);
result=test_close(result,'cell inflate',cell_inflated{2},B,1e-15,1e-15,...
                  'inflating numeric cell entries must leave them unchanged');
cell_block=blkdiag({A},{B});
result=test_close(result,'cell blkdiag upper block',cell_block{1,1},A,1e-15,1e-15,...
                  'cell blkdiag must keep the upper diagonal block');
result=test_close(result,'cell blkdiag lower block',cell_block{2,2},B,1e-15,1e-15,...
                  'cell blkdiag must keep the lower diagonal block');
if ~isempty(cell_block{1,2})||~isempty(cell_block{2,1})
    error('FAILED: cell blkdiag off-diagonal cells must be empty.');
end
result.messages{end+1}='PASS: cell blkdiag off-diagonal cells are empty.';

% Exercise recursive structure arithmetic
struct_a.alpha=[1;2];
struct_a.beta.gamma=[3 4];
struct_a.beta.delta={A,B};
struct_b.alpha=[5;6];
struct_b.beta.gamma=[7 8];
struct_b.beta.delta={C,D};
struct_c=struct_a+struct_b;
result=test_close(result,'struct plus top field',struct_c.alpha,[6;8],1e-15,1e-15,...
                  'struct plus must recurse into top-level numeric fields');
result=test_close(result,'struct plus nested field',struct_c.beta.gamma,[10 12],1e-15,1e-15,...
                  'struct plus must recurse into nested numeric fields');
result=test_close(result,'struct plus nested cell',struct_c.beta.delta{1},A+C,1e-15,1e-15,...
                  'struct plus must delegate nested cell arithmetic');
struct_d=2*struct_a;
result=test_close(result,'struct mtimes top field',struct_d.alpha,2*struct_a.alpha,1e-15,1e-15,...
                  'struct mtimes must recurse into top-level numeric fields');
result=test_close(result,'struct mtimes nested cell',struct_d.beta.delta{2},2*B,1e-15,1e-15,...
                  'struct mtimes must delegate nested cell multiplication');

% Exercise the double inflate no-op dispatch
result=test_close(result,'double inflate no-op',inflate(A),A,1e-15,1e-15,...
                  'double inflate must leave dense numeric arrays unchanged');

end


