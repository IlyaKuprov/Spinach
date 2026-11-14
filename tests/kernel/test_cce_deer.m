% Strict DEER frequency-window and species selection regression.
% Syntax: result=test_cce_deer()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_deer()

% Construct exact angular-frequency boundary values without unit conversion noise
result=new_test_result('kernel/cce_deer','CCE DEER selection',...
                       'Check strict frequency windows, isotope filters and empty selections.');
model.central_ids=1; model.bath_ids=2:6; model.multiplicities=2*ones(1,6);
model.isotopes={'E','1H','1H','1H','1H','13C'}; model.magnet_T=0;
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
fprintf('DEER exact-boundary cases passed: selected sites [3 4], [3 6], and empty.\n');

end


