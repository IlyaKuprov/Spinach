% Test: kernel/clexpn_options | Native CCE inputs and bath construction
% Deterministic lattice, occupancy, geometry, and frame tests.
%
%                         result=test_clexpn_options()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_options()

% Initialise standard harness fields only
result=new_test_result('kernel/clexpn_options','Native CCE inputs and bath construction',...
                       'Nonorthogonal cells, occupancy, density and covariance.');

% Compare nonorthogonal replication against direct Cartesian arithmetic
s=struct('mode','lattice','isotopes',{{'1H','13C'}},'magnetic',1,...
         'cell_A',[2 0 0;1 3 0;0.2 0.4 4],'fractional',[0 0 0;0.25 0.5 0.75],...
         'repeats',[2 1 1]);
[sys,~,bath]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
expected=[s.fractional*s.cell_A;(s.fractional+[1 0 0])*s.cell_A];
result=test_close(result,'Nonorthogonal row-vector cell convention',bath.coordinates_A,expected,0,0,...
                  'Nonorthogonal row-vector cell convention');
result=test_true(result,'Unique physical site IDs',numel(sys.isotopes)==4&&numel(unique(bath.ids))==4,...
                 'Unique physical site IDs');

% Check isotope occupancy without mutating the global random stream
s.species={'1H','2H'};
s.occupancy=[0.4 0.6];
s.seed=791;
state=rng;
[~,~,a]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
[~,~,b]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
result=test_true(result,'Seeded isotope occupancy is reproducible and local',isequal(a,b)&&isequal(state,rng),...
                 'Seeded isotope occupancy is reproducible and local');

% Verify active frame rotation of coordinates and Cartesian interactions
q=struct('mode','explicit','isotopes',{{'1H','13C'}},'magnetic',1,...
         'coordinates_A',[1 2 3;4 5 6],'rotation',[0 -1 0;1 0 0;0 0 1],...
         'translation_A',[2 3 4]);
q.pairs_Hz=struct('ids',[1 2],'tensor',[1 2 3;4 5 6;7 8 10]);
[~,inter,bath]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_close(result,'Coordinate covariance',bath.coordinates_A,...
                  q.coordinates_A*q.rotation.'+q.translation_A,0,0,'Coordinate covariance');
result=test_close(result,'Tensor covariance',inter.coupling.matrix{1,2},...
                  q.rotation*q.pairs_Hz.tensor*q.rotation.',0,0,'Tensor covariance');

% Preserve the preceding geometry fixture while checking tensor conversion
saved_q=q;

% Round-trip local total-frequency tensors through native Zeeman conventions
q=struct('mode','explicit','isotopes',{{'1H','E','15N'}},'magnetic',1,...
         'coordinates_A',[0 0 0;1 0 0;2 0 0]);
q.rotation=[0 -1 0;1 0 0;0 0 1];
q.tensors_Hz=cat(3,[100 2 3;2 200 4;3 4 300],diag([400 500 600]),diag([-70 90 110]));
[sys,inter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
ss=create(sys,inter);
tensor_error=0;
for k=1:3
    reference=2*pi*q.rotation*q.tensors_Hz(:,:,k)*q.rotation.';
    result=test_close(result,'Local Hz tensor native Zeeman roundtrip',...
                      ss.inter.zeeman.matrix{k},reference,1e-7,1e-10,...
                      'Nuclear ppm and electron g input must reproduce the rotated total angular-frequency tensor.');
    tensor_error=max(tensor_error,norm(ss.inter.zeeman.matrix{k}-reference,'fro'));
end
q.magnetic=0;
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch exception
    rejected=strcmp(exception.identifier,'spin_bath:zero_field');
end
result=test_true(result,'Nonzero local Hz tensor at zero field refused',rejected,...
                 'Native Zeeman inputs cannot represent a field-independent nonzero tensor at zero field.');

% Preserve representable zero tensors at zero magnetic field
q.tensors_Hz=zeros(3,3,3);
[sys,inter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
ss=create(sys,inter);
result=test_true(result,'Zero local Hz tensors at zero field accepted',...
                 all(cellfun(@(tensor) norm(tensor,'fro')==0,ss.inter.zeeman.matrix)),...
                 'The zero-field refusal is limited to nonzero tensors.');

fprintf('Local Hz roundtrip max rad/s error %.3g.\n',tensor_error);
q=saved_q;

% Select a shell through the documented centre-coordinate field
shell=struct('mode','explicit','isotopes',{{'1H','1H','1H'}},...
             'magnetic',0,'coordinates_A',[0 0 0;1 0 0;2 0 0]);
shell.selection=struct('centres_A',[0 0 0],'inner_A',0.5,'outer_A',1.5);
selected=[];
accepted=false;
try
    [~,~,selected_bath]=spin_bath(struct('magnet',shell.magnetic),struct(),rmfield(shell,'magnetic'));
    selected=selected_bath.ids;
    accepted=true;
catch

    % Convert a validator rejection into a named regression failure
end
result=test_true(result,'Documented shell centres accepted',accepted,...
                 'selection.centres_A must be accepted by the shell validator.');
result=test_true(result,'Documented shell selects physical IDs',isequal(selected,2),...
                 'The shell [0.5,1.5] around the origin must retain only site 2.');

% Test a two-dimensional number density and finite packing exclusion
r=struct('mode','random','isotopes',{{'E'}},'magnetic',0,...
         'bounds_A',[0 0 0;100 200 0],'dimension',2,'density',1e17,...
         'seed',19,'max_attempts',10000,'exclusion_A',1);
[~,~,bath]=spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
surface_density=r.density;
result=test_true(result,'Surface density gives twenty unique sites',...
                 numel(bath.ids)==20&&size(unique(bath.coordinates_A,'rows'),1)==20,...
                 'Surface density gives twenty unique sites');
result=test_close(result,'SI density and area conversion',bath.realized_density,r.density,1,1e-14,...
                  'SI density and area conversion');

% Reject duplicate positions and remove sites before preserving group closure
q=rmfield(q,{'rotation','translation_A','pairs_Hz'});
q.coordinates_A=[0 0 0;1 0 0];
q.groups={[1 2]};
q.remove_ids=2;
[~,~,bath]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_true(result,'Defect removal preserves physical group IDs',...
                 isequal(bath.ids,1)&&isequal(bath.groups,{1}),'Defect removal preserves physical group IDs');
q=rmfield(q,'remove_ids');
q.coordinates_A=zeros(2,3);
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Duplicate physical sites are refused',rejected,'Duplicate physical sites are refused');

% Check explicit physical grouping and the native isotope row contract
q=struct('mode','explicit','isotopes',{{'1H','13C'}},'magnetic',1,...
         'coordinates_A',[0 0 0;1 0 0],'groups',{{1,2}});
[sys,inter,b]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_true(result,'Complete groups and isotope row convention',...
                 isequal(b.groups,{1,2})&&isrow(sys.isotopes),...
                 'User groups explicitly partition every retained site.');

% Check scalar isotropic interactions match explicitly supplied tensors
q.pairs_Hz=struct('ids',[1 2],'tensor',37);
[~,scalar]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
q.pairs_Hz.tensor=37*eye(3);
[~,tensor]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_close(result,'Scalar isotropic tensor',scalar.coupling.matrix{1,2},...
                  tensor.coupling.matrix{1,2},0,0,'Both inputs are in Hz.');

% Check bath-to-Spinach point-dipole covariance without a Hamiltonian shortcut
q=rmfield(q,'pairs_Hz');
R=[0 -1 0;1 0 0;0 0 1];
q.rotation=R;
q.translation_A=[3 4 5];
[rsys,rinter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
rsys.output='hush';
sys.parallel={'local',1};
rsys.parallel={'local',1};
plain=create(sys,inter);
rotated=create(rsys,rinter);
T=plain.inter.coupling.matrix{1,2};
Tr=rotated.inter.coupling.matrix{1,2};
result=test_close(result,'Native dipolar transform covariance',Tr,R*T*R.',...
                  100*eps*norm(T,'fro'),100*eps,'Rigid translations leave dipoles invariant.');

% Reject invalid original site labels rather than silently dropping them
q=rmfield(q,{'rotation','translation_A'});
q.groups={[1 3]};
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Unknown group site refused',rejected,'Original candidate IDs are checked.');
q.groups={1,2};
q.pairs_Hz=struct('ids',[1 3],'tensor',eye(3));
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Unknown pair site refused',rejected,'A misspelled interaction cannot disappear.');

% Check the linear density boundary independently of the areal case
r=struct('mode','random','isotopes',{{'1H'}},'magnetic',0,'seed',12,...
         'dimension',1,'bounds_A',[0 2 3;10 2 3],'density',3e9,'max_attempts',20);
[~,~,line]=spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
result=test_true(result,'Linear density and frozen coordinates',numel(line.ids)==3&&...
                 all(line.coordinates_A(:,2)==2)&&all(line.coordinates_A(:,3)==3),...
                 'Three spins in 10 Angstrom at 3e9 per metre.');

% Expose unattainable random packing at the caller ceiling
r.count=2;
r=rmfield(r,'density');
r.exclusion_A=20;
r.max_attempts=2;
rejected=false;
try
    spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
catch e
    rejected=strcmp(e.identifier,'spin_bath:packing');
end
result=test_true(result,'Packing ceiling exposed',rejected,'No hidden relaxation of exclusion radius.');

% Validate joint assignments against count and density-derived random candidates
q=struct('mode','random','isotopes',{{'1H'}},'magnetic',1,'bounds_A',[0 0 0;1 1 1],...
         'count',3,'dimension',3,'seed',17,'max_attempts',100);
q.configurations=struct('isotopes',{{'1H','','13C'},{'','1H','1H'}},'weight',{0.3,0.7});
for density_mode=[false true]
    trial=q;
    if density_mode
        trial=rmfield(trial,'count');
        trial.density=3e30;
    end
    accepted=false;
    try
        [~,~,joint]=spin_bath(struct('magnet',trial.magnetic),struct(),rmfield(trial,'magnetic'));
        accepted=size(joint.candidate_coordinates_A,1)==3&&...
                 isequal(joint.configurations(1).ids,[1 3])&&...
                 isequal(joint.configurations(2).ids,[2 3])&&...
                 isequal([joint.configurations.weight],[0.3 0.7]);
    catch exception
        fprintf('Random joint density=%d: %s\n',density_mode,exception.message);
    end
    result=test_true(result,sprintf('Random joint candidate count density=%d',density_mode),accepted,...
                     'Each joint isotope row covers every candidate before vacancies or exclusions are applied.');

    % Refuse singleton joint rows before any candidate-mask indexing
    trial.configurations=struct('isotopes',{{'1H'}},'weight',1);
    rejected=false;
    try
        spin_bath(struct('magnet',trial.magnetic),struct(),rmfield(trial,'magnetic'));
    catch exception
        rejected=contains(exception.message,'Joint configuration size mismatch');
    end
    result=test_true(result,sprintf('Random joint singleton refused density=%d',density_mode),rejected,...
                     'Malformed joint assignments must fail validation rather than logical-mask indexing.');
end

% Preserve original candidate IDs when random sites are geometrically excluded
q.isotopes={'1H','1H','1H'};
q.remove_ids=2;
accepted=false;
try
    [~,~,joint]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
    accepted=isequal(joint.configurations(1).ids,[1 3])&&...
             isequal(joint.configurations(2).ids,3)&&sum([joint.configurations.weight])==1;
catch exception
    fprintf('Excluded random joint: %s\n',exception.message);
end
result=test_true(result,'Random joint exclusion keeps candidate IDs',accepted,...
                 'Joint rows use the pre-exclusion candidate count even when explicit isotope labels are supplied.');

% Use the same candidate count for site-resolved independent probabilities
q=rmfield(q,{'configurations','remove_ids'});
q.isotopes={'1H'};
q.species={'1H','13C'};
q.occupancy=[1 0;0 1;1 0];
accepted=false;
try
    [sampled,~,~]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
    accepted=isequal(sampled.isotopes,{'1H','13C','1H'});
catch exception
    fprintf('Random site probabilities: %s\n',exception.message);
end
result=test_true(result,'Random site-resolved occupancy count',accepted,...
                 'Random per-site independent probabilities have one row per candidate, just like joint configurations.');
fprintf('R5587 count and density each produce three candidates; retained joint IDs [1 3] and [3].\n');

fprintf('CCE bath: lattice error %.3g; density %.16g m^-2; all geometric checks complete.\n',...
        max(abs(expected-a.candidate_coordinates_A),[],'all'),surface_density);

% Validate the new ownership boundary without invoking an old public entry point
sys=struct('isotopes',{{'1H','13C'}},'magnet',1,'output','hush',...
           'disable',{{'hygiene'}},'parallel',{{'local',1}});
system=create(sys,struct());
p=struct('observed',1,'transition',struct('spins',1,'ms',[0.5 -0.5]),...
         'times',[0 1e-5 2e-5]);
[resolved,policy]=cce_options(system,p,@sin,'labframe');
result=test_true(result,'Declared physical defaults',...
                 strcmp(resolved.cce.approximation,'generalized')&&...
                 strcmp(resolved.cce.bath_state.mode,'mixed_product')&&...
                 strcmp(resolved.cce.meanfield,'none')&&policy.order_requires_convergence&&...
                 ~policy.physical_tensor_conversion&&isequal(resolved.cce.bath_spins,2),...
                 'Default resolution records assumptions and does not claim convergence.');
explicit=resolved;
explicit.cce.order=1;
explicit.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e7);
explicit.cce.approximation='projected';
explicit.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2])}});
explicit.bath_jumps=struct('support',2,'operator','L+','rate',3000);
[checked,policy]=cce_options(system,explicit,@sin,'labframe');
result=test_true(result,'Explicit physical input retained',...
                 isequal(checked.cce.bath_state.rho,explicit.cce.bath_state.rho)&&...
                 isequal(checked.bath_jumps,explicit.bath_jumps)&&...
                 isequal(checked.times,p.times)&&~policy.order_requires_convergence&&...
                 ~policy.resources_require_resolution&&~policy.accuracy_requires_resolution,...
                 'Populations and inverse-second jump rates are neither normalized nor converted.');

% Reject obsolete ownership, unknown nested fields, and invalid physical values
invalid=cell(1,12);
invalid{1}=p;
invalid{1}.observed_spins=1;
invalid{2}=p;
invalid{2}.assumptions='labframe';
invalid{3}=p;
invalid{3}.cce=struct('aproximation','projected');
invalid{4}=explicit;
invalid{4}.cce.bath_spins=1;
invalid{5}=explicit;
invalid{5}.cce.bath_state.rho={diag([1.2 -0.2])};
invalid{6}=explicit;
invalid{6}.bath_jumps.rate=-1;
invalid{7}=explicit;
invalid{7}.bath_jumps=rmfield(invalid{7}.bath_jumps,'rate');
invalid{8}=explicit;
invalid{8}.bath_jumps.support=system.comp.nspins+1;
invalid{9}=p;
invalid{9}.spins=1;
invalid{10}=p;
invalid{10}.times=[0 0];
invalid{11}=p;
invalid{11}.cce.bath_state=struct('mode','sampled','seed',NaN);
invalid{12}=p;
invalid{12}.transition.spins=2;
for k=1:numel(invalid)
    rejected=false;
    try
        cce_options(system,invalid{k},@sin,'labframe');
    catch exception
        rejected=startsWith(exception.identifier,'cce_options:');
    end
    result=test_true(result,sprintf('New-contract rejection %d',k),rejected,...
                     'The validator rejects unsupported ownership and unphysical input.');
end

% Keep frequency units explicit and retain the supplied ensemble probabilities
pump=p;
pump.pump=struct('sites',2,'center',100,'halfwidth',200,...
                 'units','Hz','fraction',0.25,'allow_empty',false);
checked=cce_options(system,pump,@sin,'labframe');
result=test_true(result,'Pump Hz retained',isequal(checked.pump,pump.pump),...
                 'No frequency conversion occurs at the option boundary.');
pump.pump.units='MHz';
rejected=false;
try
    cce_options(system,pump,@sin,'labframe');
catch
    rejected=true;
end
result=test_true(result,'Pump unit ambiguity refused',rejected,'Unsupported frequency units are not guessed.');
pump.pump.units='rad/s';
pump.pump=rmfield(pump.pump,'sites');
pump.pump.species={'13C'};
[checked,policy]=cce_options(system,pump,@sin,'labframe');
result=test_true(result,'Angular-frequency and species filters retained',...
                 isequal(checked.pump,pump.pump)&&strcmp(policy.pump_frequency_units,'rad/s'),...
                 'Explicit units and optional intersecting site/species filters belong to the experiment.');
sampled=p;
sampled.cce.bath_state=struct('mode','sampled','basis',{{eye(2)}},...
                              'probabilities',{{[0.8 0.2]}},'samples',7,'seed',19);
state_before=rng;
checked=cce_options(system,sampled,@sin,'labframe');
result=test_true(result,'Seed and probabilities retained',...
                 isequal(state_before,rng)&&isequal(checked.cce.bath_state,sampled.cce.bath_state),...
                 'Option validation does not sample or alter the global random stream.');
sampled.cce.bath_state.seed=2^32;
rejected=false;
try
    cce_options(system,sampled,@sin,'labframe');
catch
    rejected=true;
end
result=test_true(result,'Out-of-range seed refused',rejected,'The declared stream seed is exact.');

% Append sites while preserving all existing native physical input
native=struct('isotopes',{{'1H'}},'magnet',1,'labels',{{'observed'}});
inter=struct();
inter.coordinates={[0 0 0]};
inter.zeeman.scalar={3};
inter.coupling.scalar={7};
spec=struct('mode','explicit','isotopes',{{'13C'}},'coordinates_A',[1 2 3]);
[augmented,inter,bath]=spin_bath(native,inter,spec);
result=test_true(result,'Native append ownership',...
                 isequal(augmented.isotopes,{'1H','13C'})&&augmented.magnet==native.magnet&&...
                 inter.zeeman.scalar{1}==3&&inter.coupling.scalar{1,1}==7&&...
                 isequal(inter.coordinates,{[0 0 0],[1 2 3]})&&isequal(bath.spins,2)&&...
                 isequal(bath.spin_groups,{2}),...
                 'Existing tensors and labels survive; generated sites acquire physical indices.');
spec.magnetic=2;
rejected=false;
try
    spin_bath(native,inter,spec);
catch exception
    rejected=strcmp(exception.identifier,'spin_bath:ownership');
end
result=test_true(result,'Duplicate field ownership refused',rejected,'sys.magnet is authoritative.');

% Generalised density experiments supply a state rather than a transition
physical=rmfield(p,'transition');
physical.rho0=diag([0.8 0.2]);
[checked,policy]=cce_options(system,physical,@sin,'labframe');
result=test_true(result,'Explicit generalized state needs no transition',...
                 isequal(checked.rho0,physical.rho0)&&~isfield(policy,'transition'),...
                 'The experiment retains responsibility for propagation and detection.');
physical.cce.approximation='projected';
rejected=false;
try
    cce_options(system,physical,@sin,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_options:transition');
end
result=test_true(result,'Projected manifolds must be specified',rejected,...
                 'An explicit generalized density does not define a secular projection.');

% Preserve deterministic product-state mean fields as well as sampled ones
physical.cce.approximation='generalized';
physical.cce.meanfield='static';
physical.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2])}});
checked=cce_options(system,physical,@sin,'labframe');
result=test_true(result,'Static mean field with explicit product density',...
                 strcmp(checked.cce.meanfield,'static')&&isequal(checked.cce.bath_state,physical.cce.bath_state),...
                 'Static expectations need not be obtained by stochastic sampling.');
physical.cce.meanfield='none';
physical.cce.meanfield_radius=10;
rejected=false;
try
    cce_options(system,physical,@sin,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_options:meanfield');
end
result=test_true(result,'Unused mean-field radius refused',rejected,...
                 'A physical cutoff must not be silently ignored.');

end


