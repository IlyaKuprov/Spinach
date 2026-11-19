% Native ensemble ordering, elementwise density, and partitioned-cluster tests.
% Syntax:
%
%                    result=test_clexpn_ensemble()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_ensemble()

% Specify the physical system and regression targets
result=new_test_result('kernel/clexpn_ensemble','Native ensemble reduction',...
                       'Weighted estimators, elementwise density and partitioning.');
sys.magnet=0;
sys.isotopes={'E','1H','1H'};
sys.output='hush';
sys.parallel={'local',2};
inter.zeeman.scalar={0,0,0};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 2]);
inter.coupling.matrix{1,3}=diag([0 0 3]);
s=create(sys,inter);
p.observed=1;
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.times=[0 0.07 0.15];
p.cce.order=1;
p.cce.bath_state=struct('mode','sampled','basis','uniform','probabilities','uniform',...
                        'samples',8,'seed',7);
p.cce.numerics=struct('atol',1e-10,'rtol',1e-10,'max_memory_bytes',2^26);
[a,da]=clexpn(s,@free_evolution,p,'labframe');
p.cce.averaging='internal_then_expansion';
[b,db]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Noncommuting estimators',max(abs(a-b))>1e-4,...
                 'Two nontrivial state ensembles do not commute with products.');
manual=zeros(size(a));
for n=1:numel(da.reduction.samples)
    manual=manual+da.ensemble.weights(n)*da.reduction.samples{n}.value;
end
result=test_close(result,'After-expansion weighting',a,manual,1e-13,1e-13,'Complete estimates first.');
result=test_true(result,'Internal-first one reduction',isscalar(db.reduction.samples),...
                 'Average reference and cluster records before nonlinear reduction.');
p.cce.bath_state=struct('mode','mixed_product');
p.cce.order=2;
p.output='density';
[density,dd]=clexpn(s,@free_evolution,p,'labframe');
full=dd.raw{1}(end).value;
result=test_close(result,'Full-order density',density,full,1e-12,1e-12,...
                  'The exact full-bath density overrides unresolved elementwise ratios.');
for n=1:numel(p.times)
    result=test_close(result,'Density trace',trace(reshape(density(n,:),2,2)),1,1e-12,1e-12,...
                      'Raw full-bath evolution is trace preserving.');
end
p=rmfield(p,'output');
p.cce.expansion='partition';
p.cce.order=1;
p.cce.groups={[2 3]};
p.cce.partition=struct('mode','user','max_dimension',4);
[partition,dp]=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'One complete group',partition,dp.raw{1}(end).value,1e-12,1e-12,...
                  'Partition order counts groups, not physical spins.');

% Unequal configuration weights distinguish the two nonlinear estimators
records=da.raw(1:2);
records{1}(1).weight=0.25;
records{2}(1).weight=0.75;
for k=1:numel(records{1})
    records{1}(k).value(:)=k;
    records{2}(k).value(:)=2*k-1;
    records{1}(k).valid(:)=true;
    records{2}(k).valid(:)=true;
    records{1}(k).abs_error_bound(:)=0;
    records{2}(k).abs_error_bound(:)=0;
end
local=da.parameters;
local.cce.averaging='after_expansion';
after=cce_reduce(records,da.clusters,local);
local.cce.averaging='internal_then_expansion';
inside=cce_reduce(records,da.clusters,local);
result=test_true(result,'Unequal configuration weights',any(abs(after.value-inside.value)>0.1),...
                 'A nontrivial two-configuration ensemble distinguishes averaging order.');

% Truncated elementwise density keeps unresolved zero-reference masks
p=rmfield(p,'transition');
p.rho0=diag([1 0]);
p.output='density';
p.cce=rmfield(p.cce,{'groups','partition'});
p.cce.expansion='standard';
p.cce.order=1;
[~,truncated]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Truncated density masks',any(~truncated.reduction.valid,'all'),...
                 'Zero-reference elementwise ratios are reported, not repaired.');
p.cce.numerics.zero_policy='error';
caught=false;
try
    clexpn(s,@free_evolution,p,'labframe');
catch exception
    caught=strcmp(exception.identifier,'cce_reduce:unresolved');
end
result=test_true(result,'Explicit unresolved error policy',caught,...
                 'Native reduction honours the requested error policy.');

% Error bounds and masks remain part of the weighted raw-record contract
raw=da.raw;
raw{1}(2).valid(2)=false;
local=da.parameters;
masked=cce_reduce(raw,da.clusters,local);
result=test_true(result,'Mask propagation',~masked.valid(2),'An invalid weighted sample is not discarded.');
raw=da.raw;
raw{1}(2).abs_error_bound(:)=0.01;
bounded=cce_reduce(raw,da.clusters,local);
result=test_true(result,'Error propagation',any(~bounded.valid),'Reported uncertainty cannot be erased.');

% Preserve the sampled full-bath standard error and unavailable estimates
request=da.parameters;
request.cce.order=2;
request.cce.averaging='internal_then_expansion';
[~,sampled]=clexpn(s,@free_evolution,request,'labframe');
values=cellfun(@(raw)raw(end).value,sampled.raw,'UniformOutput',false);
values=cat(3,values{:});
count=size(values,3);
expected=sqrt(sum(abs(values-mean(values,3)).^2,3)/(count*(count-1)));
result=test_close(result,'Internal full-bath standard error',...
                  sampled.reduction.sampling_standard_error,expected,0,0,...
                  'Use the individual exact sampled signals before internal averaging.');
request.cce.order=1;
[~,sampled]=clexpn(s,@free_evolution,request,'labframe');
result=test_true(result,'Internal truncated uncertainty is unavailable',...
                 all(isnan(sampled.reduction.sampling_standard_error),'all'),...
                 'A nonlinear product of means has no per-sample estimator.');
request.cce.bath_state.samples=1;
[~,sampled]=clexpn(s,@free_evolution,request,'labframe');
result=test_true(result,'Single draw uncertainty is unavailable',...
                 all(isnan(sampled.reduction.sampling_standard_error),'all'),...
                 'One random sample cannot estimate its variance.');

% Limit each scheduled batch without changing the physical calculation
request=da.parameters;
request.cce.parallel=struct('mode','serial','batch_size',4,'max_tasks',1);
limited=clexpn(s,@free_evolution,request,'labframe');
result=test_close(result,'Parallel task ceiling preserves signals',limited,a,0,0,...
                  'The task limit caps batches, not the cluster family.');

% Reject an interlacing density inconsistent with its sampled basis state
request=da.parameters;
request.cce.interlaced=struct('max_states',4,'reference_state',1);
model=cce_compile(s,request,'labframe');
model.clusters=da.clusters;
ensemble=da.ensemble;
ensemble.spin_config(1).bath_rho{1}=eye(2)/2;
schedule=struct('pulse_sequence',@free_evolution,'contract',da.contract);
refused=false;
try
    cce_interlace(model,da.clusters(1),ensemble,schedule,request,da.contract.evaluator);
catch exception
    refused=strcmp(exception.identifier,'cce:spin_config')&&...
            strcmp(exception.message,['Interlaced reference bath_rho must agree with its basis '...
                                      'assignment; use explicit_product for mixed states.']);
end
result=test_true(result,'Interlaced reference density consistency',refused,...
                 'Conditional sampling must use the declared reference basis state.');

% Adapt pilot groups through the native context, then validate held-out states
inter.coupling.matrix{2,3}=0.7*eye(3);
coupled=create(sys,inter);
request=dp.parameters;
exact=clexpn(coupled,@free_evolution,request,'labframe');
request.cce.groups={2,3};
request.cce.order=1;
pilot=dp.ensemble.spin_config;
heldout=pilot;
heldout.bath_rho{1}=diag([1 0]);
request.cce.partition=struct('mode','adaptive','max_dimension',4,...
                            'patches',{{[2 3]}},'patch_atol',1e-10,'patch_rtol',1e-10,...
                            'pilot_configs',pilot,'heldout_configs',heldout,...
                            'max_rounds',2,'max_probes',100);
[adapted,adaptive]=clexpn(coupled,@free_evolution,request,'labframe');
result=test_true(result,'Native adaptive partition is held-out validated',...
                 adaptive.policy.partition.complete&&...
                 isequal(adaptive.policy.partition.groups,{[2 3]})&&...
                 ~isempty(adaptive.policy.partition.diagnostics.heldout_error),...
                 'Only the frozen pilot-selected partition reaches held-out validation.');
result=test_close(result,'Adaptive full-bath group preserves signal',adapted,exact,0,0,...
                  'Adaptive grouping retains the same exact full-bath experiment.');

end


