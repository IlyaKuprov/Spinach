% Test: kernel/cce_deer | Bath-pumped Hahn spectral selection
% Strict DEER frequency-window and species selection regression. Syntax:
%
%                    result=test_cce_deer()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_deer()

% Construct exact angular-frequency boundary values without unit conversion noise
result=new_test_result('kernel/cce_deer','CCE DEER selection',...
                       'Check strict frequency windows, isotope filters and empty selections.');
sys=struct('isotopes',{{'E','1H','1H','1H','1H','13C'}},'magnet',0,...
           'output','hush','disable',{{'hygiene'}},'parallel',{{'processes',1}});
ss=create(sys,struct());
frequencies=[-10 -9 0 10 2];
for k=1:5
    ss.inter.zeeman.matrix{k+1}(3,3)=frequencies(k);
end
p.observed=1;
p.times=[0 1];
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce=struct('approximation','generalized','order',5,'meanfield','none',...
             'numerics',struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26));
p.pump=struct('species',{{'1H'}},'center',0,'halfwidth',10,...
              'boundary','open','units','rad/s','fraction',0.5,'allow_empty',false);
model=cce_compile(ss,p,'labframe');
s=selection(model,p);
result=test_true(result,'Both exact boundaries excluded and species respected',isequal(s.selected_sites,[3 4]),'Both exact boundaries excluded and species respected');
result=test_true(result,'Boundary rule is recorded',strcmp(s.boundary_rule,'strict_abs_less_than'),'Boundary rule is recorded');
p.pump=rmfield(p.pump,'species');
p.pump.sites=[2 3 5 6];
s=selection(model,p);
result=test_true(result,'Explicit site selection intersects frequency window',isequal(s.selected_sites,[3 6]),'Explicit site selection intersects frequency window');
p.pump.halfwidth=0;
refused=false;
try
    selection(model,p);
catch exception
    refused=strcmp(exception.identifier,'bath_pumped_hahn:empty');
end
result=test_true(result,'Empty DEER window requires explicit opt-in',refused,'Empty DEER window requires explicit opt-in');
p.pump.allow_empty=true;
s=selection(model,p);
result=test_true(result,'Explicit empty-selection opt-in accepted',isempty(s.selected_sites),'Explicit empty-selection opt-in accepted');

% Check canonical rad/s values are not shifted by a second gamma B term
model.spin_system.inter.magnet=2;
model.spin_system.inter.gammas=3*ones(1,6);
p.pump.halfwidth=10;
p.pump.center=0;
s=selection(model,p);
result=test_true(result,'Canonical frequencies are not converted twice',isequal(s.selected_sites,[3 6]),'Canonical frequencies are not converted twice');

% Exercise isotope selection through canonical create/prepare metadata
sys=struct('isotopes',{{'1H','13C','1H'}},'magnet',0,'output','hush',...
           'disable',{{'hygiene'}},'parallel',{{'local',1}});
inter.zeeman.scalar={0,0,0};
ss=create(sys,inter);
q=struct('observed',1,'output','density','times',[0 0.01],...
         'rho0',ones(2)/2,'orientation',[0 0 0],...
         'transition',struct('spins',1,'ms',[0.5 -0.5]));
q.cce=struct('bath_spins',[2 3],'approximation','generalized','expansion','standard',...
             'order',2,'interaction_cutoff',Inf,'retain_raw',true,...
             'graph',struct('mode','complete','score','coupling_norm'),...
             'bath_state',struct('mode','mixed_product'),'meanfield','none',...
             'numerics',struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26));
q.pump=struct('species',{{'13C'}},'center',0,'halfwidth',1,'boundary','open','units','rad/s',...
              'fraction',0.5,'allow_empty',false);
prepared=cce_compile(ss,q,'labframe');
selected=selection(prepared,q);
result=test_true(result,'Prepared model DEER isotope selection',...
                 isequal(selected.selected_sites,2),...
                 'Species filters must read the canonical isotope table produced by create.');
filtered=pumped(ss,q);
q.pump=rmfield(q.pump,'species');
q.pump.sites=2;
by_site=pumped(ss,q);
result=test_close(result,'Public DEER species and site filters agree',...
                  filtered,by_site,0,0,...
                  'The public driver must accept documented species filters without synthetic metadata.');

% Reject a non-axial canonical Zeeman field before frequency-window selection
rotated=cce_compile(ss,q,'labframe');
rotated.zeeman{2}(3,:)=[3 0 4];
axial=q;
axial.pump.center=4;
axial.pump.halfwidth=0.1;
rejected=false;
try
    selection(rotated,axial);
catch exception
    rejected=strcmp(exception.identifier,'bath_pumped_hahn:axial');
end
result=test_true(result,'DEER rejects non-axial Zeeman selection',rejected,...
                 'The full 3*Sx+4*Sz gap is 5 rad/s, not the longitudinal 4 rad/s.');

% Validate actual public lab-frame inputs before any evaluator dispatch
rotated_system=ss;
rotated_system.inter.zeeman.matrix{2}(3,:)=[3 0 4];
rejected=false;
try
    pumped(rotated_system,axial);
catch exception
    rejected=strcmp(exception.identifier,'bath_pumped_hahn:axial');
end
result=test_true(result,'Public DEER rejects non-axial Zeeman selection',rejected,...
                 'Driver preparation must preserve the transverse field and reach the selector guard.');

% Do not apply the axial requirement to unselected physical sites
axial.pump.sites=3;
axial.pump.center=0;
axial.pump.halfwidth=1;
selected=selection(rotated,axial);
result=test_true(result,'DEER axial check respects site filters',...
                 isequal(selected.selected_sites,3),...
                 'An unrelated non-axial site must not invalidate an axial selection.');

% Reject projected dynamics that discard actual observed spin transitions
transverse=ss;
transverse=cce_test_write(transverse,[1 2],diag([2 0 0]));
projected=rmfield(q,'pump');
projected.output='signal';
projected.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
projected.rho0=[1 -1i;1i 1]/2;
for method={'cce','mecce','noise'}
    projected.cce.approximation='projected';
    sequence=@free_evolution;
    if strcmp(method{1},'mecce')
        projected.cce.basis.formalism='zeeman-liouv';
    end
    if strcmp(method{1},'noise')
        projected.cce.approximation='generalized';
        projected.output='correlation';
        projected.noise=struct('terms',struct('support',2,'operator',diag([1 -1])/2),...
                               'connected',true,'ordering','ordered','manifold','alpha');
        sequence=@noise_correlation;
    end
    rejected=false;
    try
        clexpn(transverse,sequence,projected,'labframe');
    catch exception
        rejected=ismember(exception.identifier,{'cce:population_transfer','noise_correlation:manifold'});
    end
    result=test_true(result,['Projected transition rejection: ' method{1}],rejected,...
                     'Discarding a nonzero observed transition requires generalized dynamics.');

    % Exercise the same invariant-subspace guard through interlaced recursion
    if ~strcmp(method{1},'noise')
        ip=projected;
        ip.cce.bath_state=struct('mode','interlaced','basis','uniform',...
                                 'probabilities','uniform','samples',1,'seed',19);
        ip.cce.interlaced=struct('condition_sites',{{[2 3]}},'reference_state',1,'max_states',4);
        ip.cce.averaging='after_expansion';
        rejected=false;
        try
            clexpn(transverse,sequence,ip,'labframe');
        catch exception
            rejected=ismember(exception.identifier,{'cce:population_transfer','noise_correlation:manifold'});
        end
        result=test_true(result,['Interlaced transition rejection: ' method{1}],rejected,...
                         'Conditioning must not bypass the actual-Hamiltonian invariance check.');
    end
    if isfield(projected.cce,'basis')
        projected.cce=rmfield(projected.cce,'basis');
    end
end

% Preserve explicitly allowed empty selections through every native experiment
q.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
q.pump.center=100;
q.pump.halfwidth=1;
q.pump.allow_empty=true;
for method={'cce','gcce','mecce','megcce','noise'}
    q.cce.approximation='generalized';
    q.output='signal';
    q.normalized=true;
    if ismember(method{1},{'cce','mecce'})
        q.cce.approximation='projected';
    end
    if ismember(method{1},{'mecce','megcce'})
        q.cce.basis.formalism='zeeman-liouv';
    end
    if ismember(method{1},{'gcce','megcce'})
        q.normalized=false;
    end
    sequence=@free_evolution;
    if strcmp(method{1},'noise')
        q.output='correlation';
        sequence=@noise_correlation;
        q.noise=struct('connected',false,'ordering','ordered','manifold','bath_only',...
                       'terms',struct('support',2,'operator',diag([1 -1])/2));
    end
    reference=rmfield(q,'pump');
    [expected,expected_info]=clexpn(ss,sequence,reference,'labframe');
    difference=Inf;
    try
        [empty,empty_info]=pumped(ss,q);
        difference=max(abs(empty-expected),[],'all');
        passed=isequal(empty_info.reduction.valid,expected_info.reduction.valid);
    catch
        passed=false;
    end
    result=test_true(result,['Empty DEER permitted: ' method{1}],passed,...
                     'An explicitly allowed empty selection is a no-op in every method.');
    result=test_close(result,['Empty DEER no-pulse value: ' method{1}],difference,0,...
                      1e-12,1e-12,'Removing an empty pulse must preserve the no-pulse result.');
    if isfield(q.cce,'basis')
        q.cce=rmfield(q.cce,'basis');
    end
end

% Reject obsolete controls at native option, compilation, and context boundaries
q.cce.approximation='projected';
q.output='signal';
q=rmfield(q,{'noise','pump'});
q.controls=struct('bandwidth',1);
clean=rmfield(q,'controls');
canonical=cce_compile(ss,clean,'labframe');
cluster=struct('id',1,'sites',2,'physical_ids',2);
spin_config=struct('bath_ids',[2 3],'bath_rho',{{eye(2)/2,eye(2)/2}},...
                   'state_id',1,'configuration_id',1,'weight',1);
retained_system=cce_subsystem(canonical,cluster.sites,spin_config,[]);
result=test_true(result,'Native bath configuration ordering',...
                 isequal(retained_system.physical_ids,[1 2]),...
                 'The native subsystem preserves observed-first physical spin order.');
calls={@()cce_options(ss,q,@free_evolution,'labframe'),...
       @()cce_compile(ss,q,'labframe'),@()clexpn(ss,@free_evolution,q,'labframe'),...
       @()clexpn(ss,@hahn_echo,q,'labframe'),@()clexpn(ss,@bath_pumped_hahn,q,'labframe'),...
       @()clexpn(ss,@finite_drive,q,'labframe'),@()clexpn(ss,@noise_correlation,q,'labframe')};
for n=1:numel(calls)
    rejected=false;
    try
        calls{n}();
    catch exception
        rejected=contains(exception.message,'controls');
    end
    result=test_true(result,sprintf('Unsupported controls rejected path %d',n),rejected,...
                     'An accepted but unused control container must never change physics silently.');
end

% Audit ignored nested options and inactive experiment settings
base=rmfield(q,{'controls','rho0'});
base.normalized=true;
cases=repmat({base},1,17);
sequences=repmat({@free_evolution},1,17);
cases{1}.cce.numerics.propagator='krylov';
cases{2}.cce.numerics.gpu=true;
cases{3}.cce.bath_state.temperature=300;
cases{4}.cce.bath_state.samples=2;
cases{5}.cce.graph.cutoff=1;
cases{6}.sampling=struct();
cases{7}.normalized=1;
cases{8}.spectrum=struct('omega',0);
cases{9}.bath_jumps=struct('support',2,'operator',eye(2),'rate',1);
sequences{9}=@cpmg_echo;
cases{10}.refocusing=struct('spins',1,'operator',pi*[0 1;1 0]/2,'fraction',0.5,'bandwidth',1);
sequences{10}=@hahn_echo;
cases{11}.refocusing=struct('spins',1,'operator',pi*[0 1;1 0]/2,'fraction',0.5,'duration',1);
sequences{11}=@hahn_echo;
cases{12}.cce.numerics.reduction='direct';
cases{13}.cce.approximation='generalized';
cases{13}.rho0=eye(2)/2;
cases{13}.cce.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',0.1);
cases{14}.analysis=struct('fit',false,'model','exponential');
cases{15}.cce.approximation='generalized';
cases{15}.output='correlation';
sequences{15}=@noise_correlation;
cases{15}.noise=struct('connected',false,'ordering','ordered','manifold','bath_only',...
                       'terms',struct('support',2,'operator',eye(2),'weight',2));
cases{16}=cases{15};
sequences{16}=@noise_correlation;
cases{16}.noise.terms=rmfield(cases{16}.noise.terms,'weight');
cases{16}.bath_jumps=struct('support',2,'operator',eye(2),'rate',1);
sequences{16}=@cpmg_echo;
cases{17}.bath_jumps=struct('support',2,'operator',eye(2),'rate',1,'type','local');
labels={'propagator','gpu','temperature','samples','cutoff','sampling',...
        'normalized','spectrum','bath_jumps','refocusing','refocusing','reduction','Second-order','analysis','weight','bath_jumps','type'};
for n=1:numel(cases)
    rejected=false;
    try
        clexpn(ss,sequences{n},cases{n},'labframe');
    catch exception
        rejected=contains(lower(exception.message),lower(labels{n}))||...
                 (n==4&&strcmp(exception.identifier,'cce_options:state'))||...
                 (n==1&&strcmp(exception.identifier,'clexpn:propagator'))||...
                 (n==9&&strcmp(exception.identifier,'clexpn:dissipation'))||...
                 (n==16&&strcmp(exception.identifier,'clexpn:noise'));
    end
    result=test_true(result,['Ignored option rejected: ' labels{n}],rejected,...
                     'Unsupported nested options must not silently leave a successful result unchanged.');
end

% Keep optional diagnostics separate from the natural numeric observable
retained=clexpn(ss,@free_evolution,base,'labframe');
base.cce.retain_raw=false;
suppressed=clexpn(ss,@free_evolution,base,'labframe');
result=test_true(result,'Explicit raw suppression',isnumeric(suppressed),...
                 'The natural numeric result has no compulsory raw cluster container.');
result=test_close(result,'Raw retention preserves physics',suppressed,retained,...
                  0,0,'Suppressing retained diagnostics must not change the observable.');

% Retain explicit spectral and analysis services after the native experiment
post=base;
post.cce=rmfield(post.cce,'retain_raw');
post.cce.approximation='generalized';
post.output='correlation';
post.noise=struct('connected',false,'ordering','ordered','manifold','bath_only',...
                  'terms',struct('support',2,'operator',eye(2)));
post.spectrum=struct('frequency_grid',[0 1],'convention','exp(+i*omega*t)',...
                     'window','none','weights','trapezoid','sidedness','two-sided',...
                     'time_domain','positive','stationary',true);
[actual,actual_info]=clexpn(ss,@noise_correlation,post,'labframe');
reference=cce_spectrum(post.times,actual,struct('spectrum',post.spectrum));
analysis=cce_analyse(struct('times',post.times,'value',actual,'valid',actual_info.reduction.valid),...
                     struct('analysis',struct('fit',false)));
result=test_true(result,'Driver retains requested postprocessing',...
                 isstruct(analysis)&&isfield(actual_info,'spectrum'),...
                 'Explicit analysis and spectrum requests produce their documented helper outputs.');
difference=max(abs(actual_info.spectrum.value-reference.value),[],'all');
result=test_close(result,'Driver spectrum helper parity',difference,0,...
                  0,0,'The context must invoke the same explicit Fourier contract as cce_spectrum.');

% Preserve physical mean-field settings through native experiment compilation
r=base;
r.rho0=ones(2)/2;
r.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
r.cce.meanfield='static';
r.cce.meanfield_radius=Inf;
r.times=[0 0.01];
r.cce.order=2;
for method={'cce','mecce','megcce'}
    r.cce.approximation='projected';
    if strcmp(method{1},'megcce')
        r.cce.approximation='generalized';
    end
    if ~strcmp(method{1},'cce')
        r.cce.basis.formalism='zeeman-liouv';
    end
    for control=1:3

        % Compare observed, spectrally selected, and explicit bath refocusing
        r.refocusing=struct('spins',1,'operator',pi*[0 1;1 0]/2,'fraction',0.5);
        if control>=2
            r.refocusing.spins=2;
        end
        if control==2
            r.pump=struct('center',0,'halfwidth',1,'sites',2,'boundary','open','units','rad/s',...
                          'fraction',0.5,'allow_empty',false);
        elseif isfield(r,'pump')
            r=rmfield(r,'pump');
        end
        unlimited=r;
        unlimited.cce=rmfield(unlimited.cce,'meanfield_radius');
        reference=cce_compile(ss,unlimited,'labframe');
        passed=true;
        difference=Inf;
        try
            with_radius=cce_compile(ss,r,'labframe');
            difference=double(~isequaln(with_radius.pairs,reference.pairs));
        catch
            passed=false;
        end
        name=sprintf('Static radius recursive controls %s %d',method{1},control);
        result=test_true(result,[name ' accepted'],passed,...
                         'Native compilation preserves the static-field radius without legacy recompilation.');
        result=test_close(result,[name ' parity'],difference,0,1e-12,1e-12,...
                          'An explicit infinite radius must equal the unbounded static-field default.');

    end
    if isfield(r.cce,'basis')
        r.cce=rmfield(r.cce,'basis');
    end
end

% Use a polarised exterior bath so dropping its fields is not innocuous
base=rmfield(base,intersect(fieldnames(base),{'pump','noise'}));
base.rho0=ones(2)/2;
base.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
base.cce.bath_state=struct('mode','explicit_product',...
                                 'rho',{{diag([1 0]),diag([1 0])}});
base.cce.meanfield='static';
base.cce.meanfield_radius=2;
base.cce.order=1;
cluster=struct('id',1,'sites',2);
for method={'cce','gcce','mecce','megcce'}

    % Choose each implemented evaluator without changing the missing geometry
    base.cce.approximation='generalized';
    base.output='signal';
    if ismember(method{1},{'cce','mecce'})
        base.cce.approximation='projected';
    end
    base.cce.basis.formalism='zeeman-hilb';
    if ismember(method{1},{'mecce','megcce'})
        base.cce.basis.formalism='zeeman-liouv';
    end
    model=cce_compile(ss,base,'labframe');
    base=model.parameters;
    spin_config=struct('bath_ids',[2 3],'bath_rho',{base.cce.bath_state.rho},...
                       'state_id',1,'configuration_id',1,'weight',1);
    contract=cce_capabilities(@free_evolution);
    schedule=struct('pulse_sequence',@free_evolution,'contract',contract);

    % Distinguish missing observed and exterior sites from irrelevant included sites
    for missing=[model.observed_ids 3]
        partial=model;
        partial.coordinates_A=[0 0 0;1 0 0;3 0 0];
        partial.coordinates_A(missing,:)=NaN;
        rejected=false;
        try
            contract.evaluator(partial,cluster,spin_config,schedule,base);
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates')&&...
                     contains(err.message,'meanfield_radius');
        end
        result=test_true(result,sprintf('Finite radius missing site %d %s',missing,method{1}),...
                         rejected,'Each observed and exterior site needs finite coordinates.');
    end
    spatial=model;
    spatial.coordinates_A=[0 0 0;NaN NaN NaN;3 0 0];
    accepted=true;
    try
        contract.evaluator(spatial,cluster,spin_config,schedule,base);
    catch
        accepted=false;
    end
    result=test_true(result,['Finite radius ignores included-site geometry ' method{1}],...
                     accepted,'The radius does not inspect sites already included in the cluster.');

    % Driver and direct helper calls must report the same prerequisite
    master=ismember(method{1},{'mecce','megcce'});
    calls=cell(1,5+master);
    names=cell(size(calls));
    calls(1:3)={@()clexpn(ss,@free_evolution,base,'labframe'),...
                @()contract.evaluator(model,cluster,spin_config,schedule,base),...
                @()cce_subsystem(model,cluster.sites,spin_config,[])};
    names(1:3)={'driver','cluster','controls'};
    if master
        calls{6}=@()contract.evaluator(model,cluster,spin_config,schedule,base);
        names{6}='master';
    end

    % A controlled schedule must not bypass the pre-dispatch prerequisite
    driven=schedule;
    driven.controls=true;
    calls{4}=@()contract.evaluator(model,cluster,spin_config,driven,base);
    names{4}='controlled dispatch';

    % Projected enumeration must use the same checked evaluators
    ip=base;
    ip.cce.bath_state=struct('mode','interlaced','basis','uniform',...
                             'probabilities','uniform','samples',1,'seed',7);
    ip.cce.interlaced=struct('condition_sites',{{[2 3]}},'reference_state',1,'max_states',4);
    spin_config.assignment=[1 1];
    sampled=struct('spin_config',spin_config,'weights',1,'bases',{{eye(2),eye(2)}},...
                   'probabilities',{{[0.5 0.5],[0.5 0.5]}});
    calls{5}=@()cce_interlace(model,cluster,sampled,schedule,ip,contract.evaluator);
    names{5}='interlace';
    for k=1:numel(calls)
        rejected=false;
        try
            calls{k}();
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates')&&...
                     contains(err.message,'meanfield_radius');
        end
        result=test_true(result,['Finite radius geometry ' method{1} ' ' names{k}],...
                         rejected,'A finite mean-field radius requires observed and exterior coordinates.');
    end

end

% Compile native inputs before state-dependent propagation geometry is known
base.cce.approximation='generalized';
base.output='density';
base.times=[0 0.01];
base.cce.basis.formalism='zeeman-hilb';
base.cce.meanfield_radius=2;
models={ss,ss};
models{2}.inter.coordinates={[0 0 0],[1 0 0],[3 0 0]};
for geometry=1:2
    compiled=false;
    try
        events=cce_compile(models{geometry},base,'labframe');
        compiled=isstruct(events);
    catch
    end
    result=test_true(result,sprintf('Finite-radius pre-spin_config compilation %d',geometry),...
                     compiled,'Compilation requires no state-dependent geometry; propagation does.');
end

% Complete geometry remains valid for observed, pumped, finite, and bath pulses
spatial=ss;
spatial.inter.coordinates={[0 0 0],[1 0 0],[3 0 0]};
names={'named','deer','finite','bath'};
for method={'cce','gcce','mecce','megcce'}
    for control=1:numel(names)
        if control>=3&&~strcmp(method{1},'gcce')
            continue;
        end
        p=base;
        p.cce.approximation='generalized';
        p.output='signal';
        p.normalized=true;
        if ismember(method{1},{'cce','mecce'})
            p.cce.approximation='projected';
        end
        if ismember(method{1},{'mecce','megcce'})
            p.cce.basis.formalism='zeeman-liouv';
        end
        sequence=@hahn_echo;
        if strcmp(method{1},'gcce')
            p.cce.numerics.propagator='expm';
        end
        if control==2
            p.pump=struct('center',0,'halfwidth',1,'sites',2,'boundary','open','units','rad/s',...
                          'fraction',0.5,'allow_empty',false);
            if ~strcmp(method{1},'gcce')
                p.pump.halfwidth=0;
                p.pump.allow_empty=true;
            else
                p.refocusing=struct('spins',2,'axis','z','angle',pi,'fraction',0.5);
            end
        elseif control==3
            p.cce.numerics=rmfield(p.cce.numerics,'propagator');
            sequence=@finite_drive;
            p.drive=struct('spins',1,'axis','x','amplitude',1,'detuning',0,'duration',0.005);
        elseif control==4
            p.refocusing=struct('spins',2,'axis','z','angle',pi,'fraction',0.5);
        end
        label=[method{1} ' ' names{control}];
        accepted=false;
        difference=Inf;
        try
            if control==2
                [observed,details]=pumped(spatial,p);
            else
                [observed,details]=clexpn(spatial,sequence,p,'labframe');
            end
            accepted=all(details.reduction.valid,'all');
            difference=max(abs(observed(:)-1));
        catch err
            fprintf('Finite-radius %s: %s %s\n',label,err.identifier,err.message);
        end
        result=test_true(result,['Finite-radius control dispatch ' label],accepted,...
                         'Valid geometry and a finite static-field cutoff must survive compilation.');
        result=test_close(result,['Finite-radius identity control ' label],difference,0,...
                          1e-12,0,'An uncoupled zero-field observed-spin system prepared along x retains unit signal.');
        rejected=false;
        try
            if control==2
                pumped(ss,p);
            else
                clexpn(ss,sequence,p,'labframe');
            end
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates');
        end
        result=test_true(result,['Finite-radius propagation geometry ' label],rejected,...
                         'Deferring the configuration-dependent check must preserve the physical prerequisite.');
    end
end

% A finite radius excludes an actual polarised exterior coupling
field_model=cce_compile(spatial,base,'labframe');
field_model.pairs(end+1)=struct('ids',[1 3],'tensor',zeros(3),'retained',diag([0 0 4]));
field_model.parameters.cce.meanfield_radius=2;
[near,~]=cce_subsystem(field_model,2,struct(),[]);
[h_near,~,~]=cce_generators(near,field_model.parameters,[]);
field_model.parameters.cce.meanfield_radius=4;
[far,~]=cce_subsystem(field_model,2,struct(),[]);
[h_far,R,K,field_par]=cce_generators(far,field_model.parameters,[]);
result=test_close(result,'Finite radius retains the physical exterior field',...
                  full(h_far-h_near),kron(diag([1 -1]),eye(2)),0,0,...
                  'Only the exterior spin inside the declared radius contributes its measured moment.');

% A Cartesian exterior inversion reverses the subsequent physical mean field
field_par.refocusing=struct('spins',3,'axis','x','angle',pi,'fraction',0.5);
field_par.cce.numerics.propagator='expm';
[echo,~]=hahn_echo(far,field_par,h_far,R,K);
reference=zeros(numel(field_par.times),4);
for n=1:numel(field_par.times)
    U=expm(1i*diag([1 -1])*field_par.times(n)/2)*...
      expm(-1i*diag([1 -1])*field_par.times(n)/2);
    rho=U*field_par.rho0*U';
    reference(n,:)=rho(:).';
end
result=test_close(result,'Exterior inversion changes the subsequent mean field',echo,reference,...
                  1e-13,0,'Both intervals must use the corresponding physical exterior density.');

% Invalid scalar radii cannot silently disable or empty a mean field
for radius={NaN,-1,-Inf,[1 2],1i,'2'}
    base.cce.meanfield_radius=radius{1};
    rejected=false;
    try
        cce_options(ss,base,@hahn_echo,'labframe');
    catch err
        rejected=contains(err.message,'meanfield_radius');
    end
    result=test_true(result,'Invalid finite-radius declaration',rejected,...
                     'meanfield_radius must be a nonnegative real scalar or positive Inf.');
    rejected=false;
    try
        cce_compile(ss,base,'labframe');
    catch err
        rejected=contains(err.message,'meanfield_radius');
    end
    result=test_true(result,'Invalid pre-spin_config radius declaration',rejected,...
                     'Deferring geometry must not accept an invalid cutoff during compilation.');
end

fprintf('DEER exact-boundary cases passed: selected sites [3 4], [3 6], and empty.\n');

end

% Obtain experiment-owned selection without propagation or mean-field closure
function selected=selection(model,parameters)

% Resolve the pump selection at zero propagation time
parameters.times=0;
parameters.output='signal';
parameters=rmfield(parameters,intersect(fieldnames(parameters),{'refocusing','noise'}));
parameters.cce=model.parameters.cce;
parameters.cce.meanfield='none';
parameters.cce.approximation='generalized';
if isfield(parameters.cce,'meanfield_radius')
    parameters.cce=rmfield(parameters.cce,'meanfield_radius');
end
model.parameters=parameters;
model.formalism='zeeman-hilb';
[system,~]=cce_subsystem(model,model.bath_ids,struct(),[]);
[~,selected]=bath_pumped_hahn(system,parameters,zeros(system.dimension),0,0);

end

% Apply the selected single bath pulse through the registered native experiment
function [answer,details]=pumped(spin_system,parameters)

% Map the selected bath pulse to an explicit native experiment
model=cce_compile(spin_system,parameters,'labframe');
selected=selection(model,parameters);
parameters=rmfield(parameters,'pump');
sequence=@free_evolution;
if isfield(parameters,'noise')
    sequence=@noise_correlation;
elseif ~isempty(selected.selected_sites)
    sequence=@hahn_echo;
    if ~isfield(parameters,'refocusing')
        parameters.refocusing=struct('spins',selected.selected_sites,...
                                     'axis','x','angle',pi,'fraction',0.5);
    end
    if strcmp(parameters.cce.approximation,'generalized')&&...
       (~isfield(parameters.cce,'basis')||strcmp(parameters.cce.basis.formalism,'zeeman-hilb'))
        parameters.cce.numerics.propagator='expm';
    end
end
if isempty(selected.selected_sites)&&~isfield(parameters,'noise')&&...
   strcmp(parameters.cce.approximation,'generalized')&&...
   (~isfield(parameters.cce,'basis')||strcmp(parameters.cce.basis.formalism,'zeeman-hilb'))
    sequence=@hahn_echo;
    parameters.cce.numerics.propagator='expm';
    parameters.refocusing=struct('spins',parameters.observed,'operator',zeros(2),'fraction',0.5);
end
[answer,details]=clexpn(spin_system,sequence,parameters,'labframe');

end


