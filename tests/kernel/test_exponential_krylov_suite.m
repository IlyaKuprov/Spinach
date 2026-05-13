% Tests exponential, Chebyshev, and Krylov numerical utilities. Syntax:
%
%                    result=test_exponential_krylov_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Arnoldi basis identities, Chebyshev coefficients,
% exponential drop boundary values, and Van Loan exponential-integral
% helpers against small closed-form references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_exponential_krylov_suite()

% Announce the test target
fprintf('TESTING: Exponential and Krylov numerical utilities\n');

% State the utility target of the test
result=new_test_result('kernel/exponential_krylov_suite',...
                       'Exponential and Krylov numerical utilities',...
                       'Krylov bases and exponential-integral helpers must reproduce closed-form matrix identities.');

% Check the Arnoldi orthogonality and projected recurrence identities
krylov_mat=diag([1 2 4]);
krylov_op=@(x)krylov_mat*x;
[V,H]=arnoldi(krylov_op,[1;2;3],2);
result=test_close(result,'arnoldi orthonormal basis',V'*V,eye(size(V,2)),...
                  1e-13,1e-13,...
                  'Arnoldi columns are orthonormal in the Euclidean inner product');
result=test_close(result,'arnoldi Hessenberg recurrence',...
                  krylov_mat*V(:,1:2),V*H,1e-13,1e-13,...
                  'A*V(:,1:n) equals V*H for the extended Hessenberg matrix');

% Check exact Arnoldi breakdown from an eigenvector initial condition
break_mat=[5 0;0 7];
[V_break,H_break]=arnoldi(@(x)break_mat*x,[1;0],4);
result=test_true(result,'arnoldi exact breakdown dimensions',...
                 (size(V_break,2)==1)&&isequal(size(H_break),[1 1]),...
                 'an invariant one-dimensional Krylov subspace is returned without zero padding');
result=test_close(result,'arnoldi exact breakdown value',H_break,5,...
                  1e-15,1e-15,...
                  'the projected operator on the invariant subspace is the eigenvalue');

% Check Chebyshev coefficients for a quadratic polynomial in T0,T1,T2
cheb_poly=@(x)2-3*x+4*(2*x.^2-1);
cheb_obs=cheb_coeff(cheb_poly,-1,1,8);
cheb_ref=[2 -3 4 0 0 0 0 0];
result=test_close(result,'cheb_coeff quadratic coefficients',cheb_obs,cheb_ref,...
                  1e-13,1e-13,...
                  'a degree-two Chebyshev polynomial has only its first three coefficients');

% Check the exponential drop against the boundary-value closed form
fall_obs=expdrop(5,2,0.4,5,3);
fall_time=linspace(0,0.4,5);
fall_scale=(5-2)/(1-exp(-3*0.4));
fall_ref=5-fall_scale+fall_scale*exp(-3*fall_time);
result=test_close(result,'expdrop boundary-value exponential',fall_obs,fall_ref,...
                  1e-14,1e-14,...
                  'the exponential drop satisfies both endpoint values and the requested decay rate');
result=test_true(result,'expdrop monotonic fall',all(diff(fall_obs)<0),...
                 'a positive drop rate from a larger value to a smaller value is strictly monotonic');

% Check the single exponential integral against a diagonal closed form
spin_system=local_spin_system();
left_mat=diag([1 2]);
mid_mat=[1 2;3 4]/10;
right_mat=diag([1/2 -3/2]);
int_time=0.35;
int_obs=expmint(spin_system,left_mat,mid_mat,right_mat,int_time);
int_ref=local_expmint_ref(left_mat,mid_mat,right_mat,int_time);
result=test_close(result,'expmint diagonal closed form',int_obs,int_ref,...
                  1e-12,1e-12,...
                  'each diagonal-matrix element integrates to a scalar complex exponential quotient');

% Check the zero-time shortcut in the exponential integral helper
int_zero=expmint(spin_system,left_mat,mid_mat,right_mat,0);
result=test_close(result,'expmint zero duration shortcut',int_zero,zeros(size(mid_mat)),...
                  1e-15,1e-15,...
                  'zero integration time returns the zero matrix');

% Check the nested exponential integral in the scalar nilpotent case
nested_time=0.25;
nested_obs=expmint2(spin_system,0,2,0,3,0,nested_time);
nested_ref=2*3*nested_time^2/2;
result=test_close(result,'expmint2 scalar nilpotent integral',nested_obs,nested_ref,...
                  1e-14,1e-14,...
                  'with zero generators the nested integral reduces to B*D*T^2/2');

end

% Minimal spin_system structure for propagator, clean_up, and report
function spin_system=local_spin_system()

spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.sys.output='hush';
spin_system.tols.small_matrix=1000;
spin_system.tols.prop_chop=1e-14;
spin_system.tols.dense_matrix=0.5;

end

% Closed-form integral for diagonal A and C matrices
function int_ref=local_expmint_ref(left_mat,mid_mat,right_mat,int_time)

int_ref=zeros(size(mid_mat));
left_freq=diag(left_mat);
right_freq=diag(right_mat);
for row=1:numel(left_freq)
    for col=1:numel(right_freq)
        omega=right_freq(col)-left_freq(row);
        if abs(omega)<sqrt(eps)
            int_ref(row,col)=mid_mat(row,col)*int_time;
        else
            int_ref(row,col)=mid_mat(row,col)*...
                (exp(1i*omega*int_time)-1)/(1i*omega);
        end
    end
end

end


