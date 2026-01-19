% Verification of expmint2 against numerical 
% integration.
%
% aditya.dev@weizmann.ac.il 
% ilya.kuprov@weizmann.ac.il

function expmint2_test()

% Random dimension and upper limit
n=10+randi(5); ul=randi(10)+rand();

% Generate random matrices
A=rand(n)+1i*rand(n); A=(A+A')/2;
C=rand(n)+1i*rand(n); C=(C+C')/2;
E=rand(n)+1i*rand(n); E=(E+E')/2;
B=rand(n)+1i*rand(n);
D=rand(n)+1i*rand(n);

% Bootstrap Spinach
spin_system=bootstrap();

% Call Spinach function
int_spinach=expmint2(spin_system,A,B,C,D,E,ul);

% Inner integrand and its integral
fun_inner=@(x)expm(1i*C*x)*D*expm(-1i*E*x);
int_inner=@(t)expm(-1i*C*t)*integral(fun_inner,0,t,'ArrayValued',true);

% Outer integrand and its integral
fun_outer=@(t)expm(1i*A*t)*B*int_inner(t);
int_outer=@(T)expm(-1i*A*T)*integral(fun_outer,0,T,'ArrayValued',true);

% Call Matlab integrator
int_matlab=int_outer(ul);

% Test the error norm
norm_diff=norm(int_spinach-int_matlab,'fro');
norm_comp=max([norm(int_spinach,'fro') norm(int_spinach,'fro')]);
if norm_diff/norm_comp>10*n*eps('double')
    error('double exponential integral test FAILED.');
else
    disp('double exponential integral test PASSED.');
end

end



