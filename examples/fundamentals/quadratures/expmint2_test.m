% Verification of expmint2 against numerical integration. Syntax:
%
%    expmint2_test(n,t)
%
% Parameters(optional):
%
%    n  - matrix dimension
%
%    t  - time parameter
%
% aditya.dev@weizmann.ac.il 
% ilya.kuprov@weizmann.ac.il

function expmint2_test(n,t)

% Set defaults
if nargin<1, n=randi(5); end
if nargin<2, t=randi(10)+rand(); end


% Report start
attempt_name='Verifying expmint2....';
fprintf('Running %s...\n',attempt_name);
fprintf('Size of matrices %i\n',n);
fprintf('Upper limit of Time %i\n',t);

% Generate random matrices
a=rand(n)+1i*rand(n); a=(a+a')/2;
c=rand(n)+1i*rand(n); c=(c+c')/2;
e=rand(n)+1i*rand(n); e=(e+e')/2;
b=rand(n)+1i*rand(n);
d=rand(n)+1i*rand(n);

% Compose dummy spin system
sys.magnet=14.1;
sys.isotopes={'1H'};
sys.output='hush';
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);


% Run expmint2
try
    r_fast=expmint2(spin_system,a,b,c,d,e,t);
catch me
    fprintf('expmint2 failed: %s\n',me.message);
    r_fast=[];
end

% Run numerical integration
try
    r_integral2=zeros(n);
    for p=1:n
        for q=1:n
            element_fun=@(tval,xval) arrayfun(@(tt,xx) integrand_element(tt,xx,a,b,c,d,e,t,p,q),tval,xval);
            r_integral2(p,q)=integral2(element_fun,0,t,0,@(tval) tval);
        end
    end
catch me
    fprintf('integral2 failed: %s\n',me.message);
    r_integral2=[];
end

% Report error norms
if ~isempty(r_fast)
    diff_norm=norm(r_fast-r_integral2);
    denom_norm=norm(r_integral2);
    if denom_norm==0
        if diff_norm==0
            rel_err=0;
        else
            rel_err=inf;
        end
    else
        rel_err=diff_norm/denom_norm;
    end
    fprintf('Analytical vs Numerical Difference: %e\n',diff_norm);
    fprintf('Relative Error: %e\n',rel_err);
    if rel_err<1e-8
        fprintf('SUCCESS: Verification Passed.\n');
    else
        fprintf('FAILURE: Verification Failed.\n');
    end
end

end

% Element integrand
function val=integrand_element(t,x,a,b,c,d,e,T,p,q)
exp1=expm(-1i*a*(T-t));
exp2=expm(-1i*c*(t-x));
exp3=expm(-1i*e*x);
mat=exp1*b*exp2*d*exp3;
val=mat(p,q);
end

