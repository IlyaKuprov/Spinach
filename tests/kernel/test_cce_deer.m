% Strict DEER frequency-window and species selection regression.
% Syntax: result=test_cce_deer()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_deer()

% Construct exact angular-frequency boundary values without unit conversion noise
result=new_test_result('kernel/cce_deer','CCE DEER selection',...
                       'Check strict frequency windows, isotope filters and empty selections.');
model.central_ids=1; model.bath_ids=2:6; model.multiplicities=2*ones(1,6);
model.spin_system.comp.isotopes={'E','1H','1H','1H','1H','13C'}; model.magnet_T=0;
model.spin_system.inter.gammas=zeros(1,6); model.zeeman=repmat({zeros(3)},1,6);
frequencies=[-10 -9 0 10 2];
for k=1:5, model.zeeman{k+1}(3,3)=frequencies(k); end
p.times=[0 1]; p.meanfield='none'; p.method='gcce';
p.sequence=struct('type','ideal','target','deer','placement','fraction','at',0.5,'axis','x','angle',pi);
p.deer=struct('enabled',true,'species','1H','detuning_center',0,'halfwidth',10);
s=cce_controls(model,[],[],[],p);
result=test_true(result,'Both exact boundaries excluded and species respected',isequal(s.selected_sites{1},[3 4]),'Both exact boundaries excluded and species respected');
result=test_true(result,'Boundary rule is recorded',strcmp(s.boundary_rule,'strict_abs_less_than'),'Boundary rule is recorded');
p.deer=rmfield(p.deer,'species'); p.deer.sites=[2 3 5 6]; s=cce_controls(model,[],[],[],p);
result=test_true(result,'Explicit site selection intersects frequency window',isequal(s.selected_sites{1},[3 6]),'Explicit site selection intersects frequency window');
p.deer.halfwidth=0; refused=false;
try
    cce_controls(model,[],[],[],p); catch exception, refused=strcmp(exception.identifier,'cce_controls:empty_deer'); end
result=test_true(result,'Empty DEER window requires explicit opt-in',refused,'Empty DEER window requires explicit opt-in');
p.deer.allow_empty=true; s=cce_controls(model,[],[],[],p);
result=test_true(result,'Explicit empty-selection opt-in accepted',isempty(s.selected_sites{1}),'Explicit empty-selection opt-in accepted');

% Check canonical rad/s values are not shifted by a second gamma B term
model.magnet_T=2; model.spin_system.inter.gammas=3*ones(1,6);
p.deer.halfwidth=10; p.deer.detuning_center=0; s=cce_controls(model,[],[],[],p);
result=test_true(result,'Canonical frequencies are not converted twice',isequal(s.selected_sites{1},[3 6]),'Canonical frequencies are not converted twice');

% Exercise isotope selection through canonical create/prepare metadata
sys=struct('isotopes',{{'1H','13C','1H'}},'magnet',0,'output','hush',...
           'disable',{{'hygiene'}},'parallel',{{'local',1}});
inter.zeeman.scalar={0,0,0}; ss=create(sys,inter);
q=struct('central_spins',1,'bath_spins',[2 3],'method','gcce',...
         'expansion','standard','quantity','density','times',[0 0.01],'order',2,...
         'central',struct('rho',ones(2)/2),'assumptions','labframe',...
         'orientation',[0 0 0],'interaction_cutoff',Inf,...
         'graph',struct('mode','complete','score','coupling_norm'),...
         'bath_state',struct('mode','mixed_product'),'meanfield','none',...
         'numerics',struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26));
q.sequence=struct('type','ideal','target','deer','placement','fraction',...
                  'at',0.5,'axis','x','angle',pi);
q.deer=struct('enabled',true,'species',{{'13C'}},'detuning_center',0,'halfwidth',1);
prepared=cce_prepare(ss,q); selected=cce_controls(prepared,[],[],[],q);
result=test_true(result,'Prepared model DEER isotope selection',...
                 isequal(selected.selected_sites{1},2),...
                 'Species filters must read the canonical isotope table produced by create.');
filtered=cce(ss,q); q.deer=rmfield(q.deer,'species'); q.deer.sites=2;
by_site=cce(ss,q);
result=test_close(result,'Public DEER species and site filters agree',...
                  filtered.value,by_site.value,0,0,...
                  'The public driver must accept documented species filters without synthetic metadata.');

% Preserve explicitly allowed empty selections through every driver compiler
q.central.alpha=[1;0]; q.central.beta=[0;1];
q.deer.detuning_center=100; q.deer.halfwidth=1; q.deer.allow_empty=true;
for method={'cce','gcce','mecce','megcce','noise'}
    q.method=method{1}; q.quantity='coherence'; q.jumps=struct([]);
    if strcmp(q.method,'noise')
        q.quantity='autocorrelation';
        q.noise=struct('connected',false,'ordering','ordered',...
                       'terms',struct('support',2,'operator',diag([1 -1])/2));
    end
    reference=q; reference.deer.enabled=false; reference.sequence=[]; expected=cce(ss,reference);
    passed=false; difference=Inf;
    try
        empty=cce(ss,q); difference=max(abs(empty.value-expected.value),[],'all');
        passed=isequal(empty.valid,expected.valid);
    catch
        passed=false;
    end
    result=test_true(result,['Empty DEER permitted: ' q.method],passed,...
                     'An explicitly allowed empty selection is a no-op in every method.');
    result=test_close(result,['Empty DEER no-pulse value: ' q.method],difference,0,...
                      1e-12,1e-12,'Removing an empty pulse must preserve the no-pulse result.');
end

fprintf('DEER exact-boundary cases passed: selected sites [3 4], [3 6], and empty.\n');

end


