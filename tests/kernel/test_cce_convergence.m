% Explicit independent convergence ladder and cap acceptance test.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_convergence()

% Declare finite-ladder convergence and invalid-coverage checks
result=new_test_result('kernel/cce_convergence','CCE convergence budgets',...
                       'Independent order/time ladders, resource exhaustion and held-out checks.');

% Construct a small analytic Ising bath with deterministic sampled states
sys.isotopes={'13C','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p=struct('central_spins',1,'bath_spins',[2 3],'method','cce',...
         'expansion','standard','quantity','coherence','times',[0 0.001 0.004],...
         'order',2,'central',struct('alpha',[1;0],'beta',[0;1]),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,...
         'graph',struct('mode','complete','score','coupling_norm'),...
         'meanfield','none','sequence',[],'averaging','after_expansion',...
         'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                             'samples',4,'seed',193),...
         'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));

% Compare product-factorized Ising order one against the exact order two
p.bath_state=struct('mode','mixed_product'); p.order=1;
order=struct('field','order','values',{{1,2}});
time=struct('field','times','values',{{p.times,[0 0.001 0.002 0.004]}});
p.convergence=struct('ladders',[order time],'atol',1e-10,'rtol',1e-10,'max_work',3);
converged=cce_converge(spin_system,p);
result=test_true(result,'independent ladder success',strcmp(converged.convergence,'converged'),...
                 'Ising factorization is exact and the added time mesh preserves shared points.');
result=test_true(result,'independent work accounting',converged.convergence_report.work==3&&...
                 numel(converged.convergence_report.comparisons)==2,...
                 'One baseline plus two independent refinements are recorded.');
result=test_true(result,'unchanged physics',isequal(converged.times,p.times)&&...
                 converged.input_parameters.order==1,...
                 'The reported baseline is not silently replaced by another physical calculation.');

% Exhaust the resource budget before every requested independent check
p.convergence.max_work=2; capped=cce_converge(spin_system,p);
result=test_true(result,'resource cap',strcmp(capped.convergence,'limit_reached'),...
                 'Passing one ladder cannot hide an unperformed second ladder.');
p.convergence.max_work=3; p.convergence.ladders(2).values{2}=[0 0.002 0.004]; refused=false;
try
    cce_converge(spin_system,p); catch exception, refused=strcmp(exception.identifier,'cce_converge:coverage'); end
result=test_true(result,'missing times refusal',refused,...
                 'No interpolation or discarded requested time may manufacture convergence.');

% Require independent held-out samples rather than a common-seed pilot alone
p.convergence.ladders=struct('field','bath_state.samples','values',{{4,8}});
p.bath_state=struct('mode','sampled','basis','uniform','probabilities','uniform','samples',4,'seed',193);
refused=false;
try
    cce_converge(spin_system,p); catch exception, refused=strcmp(exception.identifier,'cce_converge:sampling'); end
result=test_true(result,'held-out sampling requirement',refused,...
                 'A common-random-number ladder is not an independent uncertainty check.');
p.convergence.checks={struct('bath_state',p.bath_state)};
p.convergence.checks{1}.bath_state.seed=194;
p.convergence.atol=2; p.convergence.rtol=0;
held_out=cce_converge(spin_system,p);
result=test_true(result,'held-out check execution',held_out.convergence_report.work==3&&...
                 strcmp(held_out.convergence_report.comparisons(end).field,'independent_check'),...
                 'The explicit independent seed check consumes its own work unit.');

% Refuse numerical budgets too strict for a purported converged output
p.bath_state=struct('mode','mixed_product'); p.numerics.atol=realmin; p.numerics.rtol=0;
p.convergence=struct('ladders',order,'atol',1e-8,'rtol',0,'max_work',2);
invalid=cce_converge(spin_system,p);
result=test_true(result,'invalid coverage cannot converge',~strcmp(invalid.convergence,'converged'),...
                 'Identical approximations with failed numerical validity are not convergence.');
fprintf('CCE independent ladder discrepancy %.3g; work %d; capped status %s; invalid status %s.\n',...
        max([converged.convergence_report.comparisons.difference]),...
        converged.convergence_report.work,capped.convergence,invalid.convergence);

end


