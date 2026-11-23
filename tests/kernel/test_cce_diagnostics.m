% Test: kernel/cce_diagnostics | CCE numerical diagnostics
% Coverage and unmodified density physicality diagnostics acceptance test. Syntax:
%
%                    result=test_cce_diagnostics()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_diagnostics()

% Declare numerical coverage and density physicality acceptance checks
result=new_test_result('kernel/cce_diagnostics','CCE numerical diagnostics',...
                       'Coverage, error budgets, unphysical density and raw residuals.');

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
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^26;
spin_system=create(sys,inter);
p=struct('observed',1,'times',[0 0.001 0.004],...
         'transition',struct('spins',1,'alpha',[1;0],'beta',[0;1]),'orientation',[0 0 0]);
p.cce=struct('bath_spins',[2 3],'approximation','generalized',...
             'expansion','standard','order',2,'interaction_cutoff',Inf,...
             'meanfield','none','averaging','after_expansion');

% Check physical full-order density against an independent analytic matrix
p.rho0=ones(2)/2;
p.output='density';
[physical,physical_info]=clexpn(spin_system,@free_evolution,p,'labframe');
model_p=cce_options(spin_system,p,@cce_compile,'labframe');
model=cce_compile(spin_system,model_p,'labframe');
density=reshape(physical.',2,2,numel(p.times));
physical_checks=cce_diagnostics(model,physical_info.raw,physical_info.reduction,p);
reference=cos(17*pi*p.times).*cos(29*pi*p.times);
expected=zeros(2,2,numel(p.times));
for k=1:numel(p.times)
    expected(:,:,k)=[1 reference(k);reference(k) 1]/2;
end
result=test_close(result,'analytic reduced density',density,expected,1e-12,1e-12,...
                  'Independent diagonal Ising evolution gives the bath-averaged density.');
result=test_true(result,'physical density accepted',physical_checks.acceptable&&...
                 physical_checks.physicality_checked,...
                 'Complete exact density is normalized, Hermitian and positive.');

% Audit common result fields across standard, partition, ME, and noise manifolds
for method={'cce','gcce','mecce','megcce','noise'}
    q=p;
    q.output='signal';
    q.normalized=ismember(method{1},{'cce','mecce'});
    sequence=@free_evolution;
    q.cce.approximation='projected';
    variant=spin_system;
    if ismember(method{1},{'gcce','megcce','noise'})
        q.cce.approximation='generalized';
    end
    if ismember(method{1},{'mecce','megcce'})
        variant.bas.formalism='zeeman-liouv';
    end
    if strcmp(q.cce.approximation,'projected')
        q=rmfield(q,'rho0');
    end
    if strcmp(method{1},'noise')
        q.output='correlation';
        sequence=@noise_correlation;
        q.noise=struct('connected',false,'ordering','ordered','manifold','bath_only',...
                       'terms',struct('support',2,'operator',diag([1 -1])/2));
    end
    for expansion={'standard','partition'}
        q.cce.expansion=expansion{1};
        q.cce.groups={2,3};
        q.cce.partition=struct('mode','user','max_dimension',4);
        [checked,checked_info]=clexpn(variant,sequence,q,'labframe');
        checked_checks=cce_diagnostics(model,checked_info.raw,checked_info.reduction,q);
        result=test_true(result,['Common schema: ' method{1} ' ' expansion{1}],...
                         isnumeric(checked)&&isnumeric(physical)&&size(checked,1)==size(physical,1),...
                         'Every native experiment returns its natural numeric time-by-channel answer.');
        result=test_true(result,['Common diagnostics: ' method{1} ' ' expansion{1}],...
                         isequal(sort(fieldnames(checked_checks)),sort(fieldnames(physical_checks))),...
                         'Unavailable density diagnostics must be empty fields, not a different schema.');
        result=test_close(result,['Final prefix: ' method{1} ' ' expansion{1}],...
                          checked_info.reduction.by_order(:,:,end),checked,0,0,...
                          'The final by-order estimate must equal the returned estimate.');
        result=test_true(result,['Final reasons: ' method{1} ' ' expansion{1}],...
                         isequal(checked_info.reduction.by_order_reason(:,:,end),checked_info.reduction.reason),...
                         'A full-bath bypass must update the final prefix reason codes too.');
    end
end

% Retain all sampled prefix means rather than the first sampled trajectory
q=p;
q.bath_samples=4;
q.bath_seed=193;
[sampled,sampled_info]=clexpn(spin_system,@free_evolution,q,'labframe');
prefixes=zeros(size(sampled_info.reduction.by_order));
for k=1:numel(sampled_info.raw)
    records=sampled_info.raw{k};
    [records.weight]=deal(1);
    reduced=cce_reduce(model,{records},sampled_info.clusters,sampled_info.parameters);
    reduced.by_order(:,:,end)=sampled_info.raw{k}(end).value;
    prefixes=prefixes+sampled_info.ensemble.weights(k)*reduced.by_order;
end
result=test_close(result,'Sampled by-order means',sampled_info.reduction.by_order,prefixes,1e-12,1e-12,...
                  'Each displayed prefix must use the same sample averaging as the final result.');
result=test_close(result,'Sampled final prefix',sampled_info.reduction.by_order(:,:,end),sampled,...
                  1e-12,1e-12,'The result must not expose the first realization as an averaged prefix.');

% Unavailable residuals and sampling errors remain unavailable in diagnostic aggregates
records=physical_info.raw;
records{1}(1).eigen_residual=NaN;
records{1}(1).orthogonality_residual=NaN;
modified=physical_info.reduction;
modified.sampling_standard_error(:)=NaN;
unavailable=cce_diagnostics(model,records,modified,p);
result=test_true(result,'Unavailable diagnostic aggregates',...
                 isnan(unavailable.maximum_eigen_residual)&&...
                 isnan(unavailable.maximum_orthogonality_residual)&&...
                 all(isnan(unavailable.sampling_standard_error),'all'),...
                 'Missing estimates are not silently discarded by maximum or copy operations.');

% Insert a negative eigenvalue without allowing automatic matrix repair
modified=physical_info.reduction;
modified.value(2,:)=[1.2 0 0 -0.2];
before=reshape(modified.value.',2,2,numel(p.times));
diagnosed=cce_diagnostics(model,physical_info.raw,modified,p);
result=test_true(result,'negative density rejected',~diagnosed.physicality_pass&&...
                 abs(diagnosed.minimum_eigenvalue(2)+0.2)<eps,...
                 'Negative populations are reported, not clipped.');
result=test_true(result,'density untouched',isequal(reshape(modified.value.',2,2,numel(p.times)),before),...
                 'Diagnostics do not modify the reconstructed matrix.');

% Diagnose trace and Hermiticity errors independently
modified=physical_info.reduction;
modified.value(2,1)=modified.value(2,1)+0.1;
modified.value(2,3)=modified.value(2,3)+0.2i;
diagnosed=cce_diagnostics(model,physical_info.raw,modified,p);
result=test_true(result,'trace and Hermiticity',diagnosed.trace_error(2)>0.09&&...
                 diagnosed.hermiticity_error(2)>0.2&&~diagnosed.acceptable,...
                 'Independent physicality defects are preserved in diagnostics.');

% Never censor one failed requested time or a nonfinite error bound
modified=physical_info.reduction;
modified.valid(2,1)=false;
diagnosed=cce_diagnostics(model,physical_info.raw,modified,p);
result=test_true(result,'invalid time coverage',diagnosed.valid_time_coverage==2/3&&...
                 ~diagnosed.acceptable&&contains(diagnosed.invalid_reasons{2,1},'reduction_invalid'),...
                 'Every requested time participates in the coverage denominator.');
modified=physical_info.reduction;
modified.abs_error_bound(1,1)=NaN;
diagnosed=cce_diagnostics(model,physical_info.raw,modified,p);
result=test_true(result,'nonfinite bound',~diagnosed.acceptable&&...
                 contains(diagnosed.invalid_reasons{1,1},'invalid_error_bound'),...
                 'A finite value with a nonfinite error estimate is still invalid.');
modified=physical_info.reduction;
modified.abs_error_bound(1,1)=1;
diagnosed=cce_diagnostics(model,physical_info.raw,modified,p);
result=test_true(result,'numerical budget',~diagnosed.acceptable&&...
                 ~diagnosed.numerical_budget_met(1,1),...
                 'An exceeded propagation/reduction budget is not a truncation estimate.');
result=test_true(result,'raw diagnostics',physical_checks.raw_records==4&&...
                 physical_checks.raw_invalid==0&&...
                 isempty(physical_checks.truncation_error),...
                 'Raw task count and residuals remain separate from unknown truncation error.');
fprintf('CCE density error %.3g; trace %.3g; Hermiticity %.3g; minimum eigenvalue %.3g.\n',...
        max(abs(density-expected),[],'all'),max(physical_checks.trace_error),...
        max(physical_checks.hermiticity_error),min(physical_checks.minimum_eigenvalue));

end


