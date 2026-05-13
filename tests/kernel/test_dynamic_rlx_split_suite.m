% Tests relaxation-superoperator component splitting. Syntax:
%
%              result=test_dynamic_rlx_split_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that rlx_split() separates single-spin longitudinal,
% single-spin transverse, and multi-spin relaxation blocks without overlap.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_rlx_split_suite()

% Announce the test target
fprintf('TESTING: Relaxation component splitting\n');

% State the relaxation-splitting target of the test
result=new_test_result('kernel/dynamic_rlx_split_suite',...
                       'Relaxation component splitting',...
                       'rlx_split() must partition longitudinal, transverse, and mixed relaxation components.');

% Build a two-spin spherical-tensor Liouville basis
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.temperature=300;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Interpret the basis using the same documented state categories
[L,M]=lin2lm(spin_system.bas.basis);
sso_mask=(sum(logical(spin_system.bas.basis),2)==1);
mso_mask=(sum(logical(spin_system.bas.basis),2)>1);
long_sso_mask=any((L>0)&(M==0),2)&sso_mask;
tran_sso_mask=any((L>0)&(M~=0),2)&sso_mask;

% Build a diagonal relaxation matrix with a zero unit-state element
matrix_dim=size(spin_system.bas.basis,1);
diag_vals=(1:matrix_dim)';
diag_vals(~(long_sso_mask|tran_sso_mask|mso_mask))=0;
R=spdiags(diag_vals,0,matrix_dim,matrix_dim);

% Build the mathematically expected non-overlapping blocks
R1_ref=R;
R1_ref(~long_sso_mask,~long_sso_mask)=0;
R2_ref=R;
R2_ref(~tran_sso_mask,~tran_sso_mask)=0;
Rm_ref=R;
Rm_ref(~mso_mask,~mso_mask)=0;

% Split the relaxation superoperator with the production helper
[R1,R2,Rm]=rlx_split(spin_system,R);

% Check the three blocks against their category masks
result=test_close(result,'longitudinal block',R1,R1_ref,0,0,...
                  'R1 must contain only single-spin longitudinal states');
result=test_close(result,'transverse block',R2,R2_ref,0,0,...
                  'R2 must contain only single-spin transverse states');
result=test_close(result,'mixed block',Rm,Rm_ref,0,0,...
                  'Rm must contain only multi-spin, or otherwise mixed, states');
result=test_close(result,'block partition',R1+R2+Rm,R,0,0,...
                  'the split relaxation blocks must exactly reconstruct a diagonal category-partitioned matrix');

end


