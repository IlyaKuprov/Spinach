% Explicit independent convergence ladder and cap acceptance test. Syntax:
%
%                    result=test_cce_convergence()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_convergence()

% Declare finite-ladder convergence and invalid-coverage checks
result=new_test_result('kernel/cce_convergence','CCE convergence budgets',...
                       'Independent order/time ladders, resource exhaustion and held-out checks.');

% Construct a small analytic Ising bath with deterministic sampled states
sys.isotopes={'13C','13C','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.001 0.004];
p.orientation=[0 0 0];
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce=struct('bath_spins',[2 3],'approximation','projected','expansion','standard',...
             'order',2,'interaction_cutoff',Inf,...
             'graph',struct('mode','complete','score','coupling_norm'),...
             'meanfield','none','averaging','after_expansion',...
             'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                                 'samples',4,'seed',193),...
             'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));

% Compare product-factorised Ising order one against the exact order two
p.cce.bath_state=struct('mode','mixed_product');
p.cce.order=1;
order=struct('field','cce.order','values',{{1,2}});
time=struct('field','times','values',{{p.times,[0 0.001 0.002 0.004]}});
p.cce.convergence=struct('ladders',[order time],'atol',1e-10,'rtol',1e-10,'max_work',3);
[~,converged_details]=cce_converge(spin_system,p,@free_evolution,'labframe');
result=test_true(result,'independent ladder success',strcmp(converged_details.convergence,'converged'),...
                 'Ising factorization is exact and the added time mesh preserves shared points.');
result=test_true(result,'independent work accounting',converged_details.convergence_report.work==3&&...
                 numel(converged_details.convergence_report.comparisons)==2,...
                 'One baseline plus two independent refinements are recorded.');
result=test_true(result,'unchanged physics',isequal(converged_details.parameters.times,p.times)&&...
                 converged_details.input_parameters.cce.order==1,...
                 'The reported baseline is not silently replaced by another physical calculation.');

% Do not certify an inaccurate baseline when later ladder steps agree
inter.coupling.matrix{2,3}=diag([9 9 0]);
coupled=create(sys,inter);
q=p;
q.cce.order=1;
q.times=[0 0.01 0.04];
q.cce.bath_state=struct('mode','mixed_product');
q.cce.convergence.ladders=struct('field','cce.order','values',{{1,2,3}});
q.cce.convergence.max_work=3;
[ladder,ladder_details]=cce_converge(coupled,q,@free_evolution,'labframe');
reference=q;
reference.cce=rmfield(reference.cce,'convergence');
reference.cce.order=2;
exact=clexpn(coupled,@free_evolution,reference,'labframe');
baseline_error=max(abs(ladder-exact));
result=test_true(result,'Returned baseline must converge',...
                 strcmp(ladder_details.convergence,'not_converged')&&baseline_error>1e-3,...
                 'Agreement of orders two and three cannot certify the returned order-one baseline.');
fprintf('R1227 baseline error %.16g; ladder status %s.\n',baseline_error,ladder_details.convergence);

% Exhaust the resource budget before every requested independent check
p.cce.convergence.max_work=2;
[~,capped_details]=cce_converge(spin_system,p,@free_evolution,'labframe');
result=test_true(result,'resource cap',strcmp(capped_details.convergence,'limit_reached'),...
                 'Passing one ladder cannot hide an unperformed second ladder.');
p.cce.convergence.max_work=3;
p.cce.convergence.ladders(2).values{2}=[0 0.002 0.004];
refused=false;
try
    cce_converge(spin_system,p,@free_evolution,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_converge:coverage');
end
result=test_true(result,'missing times refusal',refused,...
                 'No interpolation or discarded requested time may manufacture convergence.');

% Require independent held-out samples rather than a common-seed pilot alone
p.cce.convergence.ladders=struct('field','cce.bath_state.samples','values',{{4,8}});
p.cce.bath_state=struct('mode','sampled','basis','uniform','probabilities','uniform','samples',4,'seed',193);
refused=false;
try
    cce_converge(spin_system,p,@free_evolution,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_converge:sampling');
end
result=test_true(result,'held-out sampling requirement',refused,...
                 'A common-random-number ladder is not an independent uncertainty check.');
p.cce.convergence.checks={struct('cce',struct('bath_state',p.cce.bath_state))};
p.cce.convergence.checks{1}.cce.bath_state.seed=194;
p.cce.convergence.atol=2;
p.cce.convergence.rtol=0;
[~,held_out_details]=cce_converge(spin_system,p,@free_evolution,'labframe');
result=test_true(result,'held-out check execution',held_out_details.convergence_report.work==3&&...
                 strcmp(held_out_details.convergence_report.comparisons(end).field,'independent_check'),...
                 'The explicit independent seed check consumes its own work unit.');

% Refuse numerical budgets too strict for a purported converged output
p.cce.bath_state=struct('mode','mixed_product');
p.cce.numerics.atol=realmin;
p.cce.numerics.rtol=0;
p.cce.convergence=struct('ladders',order,'atol',1e-8,'rtol',0,'max_work',2);
[~,invalid_details]=cce_converge(spin_system,p,@free_evolution,'labframe');
result=test_true(result,'invalid coverage cannot converge',~strcmp(invalid_details.convergence,'converged'),...
                 'Identical approximations with failed numerical validity are not convergence.');

% A partially unavailable sampling estimate cannot become a finite aggregate
request=converged_details.input_parameters;
request.cce=rmfield(request.cce,'convergence');
request.normalized=true;
request.times=[0 0.001 0.01];
request.cce.bath_state=struct('mode','sampled','samples',4,'seed',193,...
                              'basis','uniform','probabilities','uniform');
[~,available]=clexpn(spin_system,@free_evolution,request,'labframe');
request.cce.numerics.atol=mean(available.reduction.abs_error_bound([1 end]));
request.cce.numerics.rtol=0;
request.cce.convergence=struct('ladders',order,'atol',1,'rtol',0,'max_work',2);
[~,unavailable]=cce_converge(spin_system,request,@free_evolution,'labframe');
errors=unavailable.reduction.sampling_standard_error;
result=test_true(result,'Unavailable convergence sampling aggregate',...
                 any(isfinite(errors))&&any(isnan(errors))&&...
                 isnan(unavailable.convergence_report.comparisons.sampling_error),...
                 'A missing estimate at one requested time is not discarded in favour of finite times.');

% Reject malformed orders before saturating a legitimate finite-bath refinement
for value={NaN,Inf,1.5,'2'}
    request=converged_details.input_parameters;
    request.cce.convergence.ladders=struct('field','cce.order','values',{{1,value{1}}});
    invalid_order=false;
    try
        cce_converge(spin_system,request,@free_evolution,'labframe');
    catch exception
        invalid_order=strcmp(exception.identifier,'cce_converge:input');
    end
    result=test_true(result,'Malformed order ladder rejected',invalid_order,...
                     'Nonfinite, fractional, or nonnumeric orders cannot become a valid full-bath request.');
end

fprintf('CCE independent ladder discrepancy %.3g; work %d; capped status %s; invalid status %s.\n',...
        max([converged_details.convergence_report.comparisons.difference]),...
        converged_details.convergence_report.work,capped_details.convergence,invalid_details.convergence);

% Resolve the default bath before capping order refinements at its finite size
default_bath=converged_details.input_parameters;
default_bath.cce=rmfield(default_bath.cce,'bath_spins');
default_bath.cce.convergence.ladders=struct('field','cce.order','values',{{1,2,5}});
[default_value,default_details]=cce_converge(spin_system,default_bath,@free_evolution,'labframe');
explicit_bath=default_bath;
explicit_bath.cce.bath_spins=[2 3];
[explicit_value,explicit_details]=cce_converge(spin_system,explicit_bath,@free_evolution,'labframe');
result=test_close(result,'Default bath order ladder',default_value,explicit_value,0,0,...
                  'Default and explicit physical bath IDs preserve the same baseline.');
result=test_true(result,'Default bath finite-order cap',...
                 isequaln(default_details.convergence_report,explicit_details.convergence_report)&&...
                 ~isfield(default_details.input_parameters.cce,'bath_spins'),...
                 'Orders beyond the bath saturate without changing the recorded caller input.');

end


