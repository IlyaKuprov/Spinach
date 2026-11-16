% Native ensemble ordering, elementwise density and partitioned-cluster tests.
% Syntax: result=test_clexpn_ensemble()
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_ensemble()
result=new_test_result('kernel/clexpn_ensemble','Native ensemble reduction',...
                       'Weighted estimators, elementwise density and partitioning.');
sys.magnet=0; sys.isotopes={'E','1H','1H'}; sys.output='hush';
sys.parallel={'local',2}; inter.zeeman.scalar={0,0,0};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 2]);
inter.coupling.matrix{1,3}=diag([0 0 3]); s=create(sys,inter);
p.observed=1; p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.times=[0 0.07 0.15]; p.cce.order=1;
p.cce.bath_state=struct('mode','sampled','basis','uniform','probabilities','uniform',...
                        'samples',8,'seed',7);
p.cce.numerics=struct('atol',1e-10,'rtol',1e-10,'max_memory_bytes',2^26);
[a,da]=clexpn(s,@free_evolution,p,'labframe');
p.cce.averaging='internal_then_expansion'; [b,db]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Noncommuting estimators',max(abs(a-b))>1e-4,...
                 'Two nontrivial state ensembles do not commute with products.');
manual=zeros(size(a));
for n=1:numel(da.reduction.samples)
    manual=manual+da.ensemble.weights(n)*da.reduction.samples{n}.value;
end
result=test_close(result,'After-expansion weighting',a,manual,1e-13,1e-13,'Complete estimates first.');
result=test_true(result,'Internal-first one reduction',isscalar(db.reduction.samples),...
                 'Average reference and cluster records before nonlinear reduction.');
p.cce.bath_state=struct('mode','mixed_product'); p.cce.order=2;
p.output='density'; [density,dd]=clexpn(s,@free_evolution,p,'labframe');
full=dd.raw{1}(end).value;
result=test_close(result,'Full-order density',density,full,1e-12,1e-12,...
                  'The exact full-bath density overrides unresolved elementwise ratios.');
for n=1:numel(p.times)
    result=test_close(result,'Density trace',trace(reshape(density(n,:),2,2)),1,1e-12,1e-12,...
                      'Raw full-bath evolution is trace preserving.');
end
p=rmfield(p,'output'); p.cce.expansion='partition'; p.cce.order=1;
p.cce.groups={[2 3]}; p.cce.partition=struct('mode','user','max_dimension',4);
[partition,dp]=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'One complete group',partition,dp.raw{1}(end).value,1e-12,1e-12,...
                  'Partition order counts groups, not physical spins.');

% Unequal configuration weights distinguish the two nonlinear estimators
records=da.raw(1:2); records{1}(1).weight=0.25; records{2}(1).weight=0.75;
for k=1:numel(records{1})
    records{1}(k).value(:)=k; records{2}(k).value(:)=2*k-1;
    records{1}(k).valid(:)=true; records{2}(k).valid(:)=true;
    records{1}(k).abs_error_bound(:)=0; records{2}(k).abs_error_bound(:)=0;
end
local=da.parameters; local.cce.averaging='after_expansion';
after=cce_reduce(records,da.clusters,local);
local.cce.averaging='internal_then_expansion'; inside=cce_reduce(records,da.clusters,local);
result=test_true(result,'Unequal configuration weights',any(abs(after.value-inside.value)>0.1),...
                 'A nontrivial two-configuration ensemble distinguishes averaging order.');

% Truncated elementwise density keeps unresolved zero-reference masks
p=rmfield(p,'transition'); p.rho0=diag([1 0]); p.output='density';
p.cce=rmfield(p.cce,{'groups','partition'}); p.cce.expansion='standard'; p.cce.order=1;
[~,truncated]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Truncated density masks',any(~truncated.reduction.valid,'all'),...
                 'Zero-reference elementwise ratios are reported, not repaired.');
p.cce.numerics.zero_policy='error'; caught=false;
try
    clexpn(s,@free_evolution,p,'labframe');
catch exception
    caught=strcmp(exception.identifier,'cce_reduce:unresolved');
end
result=test_true(result,'Explicit unresolved error policy',caught,...
                 'Native reduction honours the requested error policy.');

% Error bounds and masks remain part of the weighted raw-record contract
raw=da.raw; raw{1}(2).valid(2)=false; local=da.parameters;
masked=cce_reduce(raw,da.clusters,local);
result=test_true(result,'Mask propagation',~masked.valid(2),'An invalid weighted sample is not discarded.');
raw=da.raw; raw{1}(2).abs_error_bound(:)=0.01;
bounded=cce_reduce(raw,da.clusters,local);
result=test_true(result,'Error propagation',any(~bounded.valid),'Reported uncertainty cannot be erased.');
end


