% Verification script for expmint2: Tests the formula against numerical
% integration

attempt_name = 'Verifying expmint2';
fprintf('Running %s...\n', attempt_name);

% Size of matrices (all square of same size)
N = randi(5);
T = randi(10); % time parameter
fprintf('Size of matrices %i \n', N);
fprintf('Upper limit of Time %i \n', T);
% 2. Generate random matrices A, C, E are Hermitian
A = rand(N) + 1i*rand(N); A = (A + A')/2;
C = rand(N) + 1i*rand(N); C = (C + C')/2;
E = rand(N) + 1i*rand(N); E = (E + E')/2;

% B, D are arbitrary
B = rand(N) + 1i*rand(N);
D = rand(N) + 1i*rand(N);

% 3. Run expmint2 spin_system dummy
sys.magnet   = 14.1;
sys.isotopes = {'1H'};
sys.output = 'hush'; % To silence report
% Chemical shift (ppm): on-resonance reference
inter.zeeman.scalar = {0.0};

% Basis set
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);



try
    R_fast = expmint2(spin_system, A, B, C, D, E, T);
catch ME
    fprintf('expmint2 failed: %s\n', ME.message);
    R_fast = [];
end

function val = integrand_element(t, x, A, B, C, D, E, T, p, q)
exp1 = expm(-1i*A*(T - t));
exp2 = expm(-1i*C*(t - x));
exp3 = expm(-1i*E*x);
mat = exp1*B*exp2*D*exp3;
val = mat(p, q);
end


try
    R_integral2 = zeros(N);
    for p = 1:N
        for q = 1:N
            element_fun = @(t, x) arrayfun(@(tt, xx) integrand_element(tt, xx, A, B, C, D, E, T, p, q), t, x);
            R_integral2(p, q) = integral2(element_fun, 0, T, 0, @(t) t);
        end
    end

catch ME
    fprintf('integral2 failed: %s\n', ME.message);
    R_integral2 = [];
end

if ~isempty(R_fast)
    diff_norm = norm(R_fast - R_integral2);
    denom_norm = norm(R_integral2);
    if denom_norm == 0
        if diff_norm == 0
            rel_err = 0;
        else
            rel_err = inf;
        end
    else
        rel_err = diff_norm / denom_norm;
    end

    fprintf('Analytical vs Numerical Difference: %e\n', diff_norm);
    fprintf('Relative Error: %e\n', rel_err);

    if rel_err < 1e-6 % Relaxed tolerance for Riemann sum
        fprintf('SUCCESS: Verification Passed.\n');
    else
        fprintf('FAILURE: Verification Failed.\n');
    end
end
