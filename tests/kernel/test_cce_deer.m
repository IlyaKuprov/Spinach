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

% Reject a non-axial canonical Zeeman field before frequency-window selection
rotated=cce_prepare(ss,q); rotated.zeeman{2}(3,:)=[3 0 4];
axial=q; axial.deer.detuning_center=4; axial.deer.halfwidth=0.1;
rejected=false;
try
    cce_controls(rotated,[],[],[],axial);
catch exception
    rejected=strcmp(exception.identifier,'cce_controls:deer_axial');
end
result=test_true(result,'DEER rejects non-axial Zeeman selection',rejected,...
                 'The full 3*Sx+4*Sz gap is 5 rad/s, not the longitudinal 4 rad/s.');

% Validate actual public lab-frame inputs before any evaluator dispatch
rotated_system=ss; rotated_system.inter.zeeman.matrix{2}(3,:)=[3 0 4];
rejected=false;
try
    cce(rotated_system,axial);
catch exception
    rejected=strcmp(exception.identifier,'cce_controls:deer_axial');
end
result=test_true(result,'Public DEER rejects non-axial Zeeman selection',rejected,...
                 'Driver preparation must preserve the transverse field and reach the selector guard.');

% Do not apply the axial requirement to unselected physical sites
axial.deer.sites=3; axial.deer.detuning_center=0; axial.deer.halfwidth=1;
selected=cce_controls(rotated,[],[],[],axial);
result=test_true(result,'DEER axial check respects site filters',...
                 isequal(selected.selected_sites{1},3),...
                 'An unrelated non-axial site must not invalidate an axial selection.');

% Reject conditional dynamics that discard actual central spin transitions
transverse=ss; transverse.inter.coupling.matrix{1,2}=diag([2 0 0]);
conditional=q; conditional=rmfield(conditional,'deer');
conditional.quantity='coherence'; conditional.sequence=[]; conditional.central.alpha=[1;0]; conditional.central.beta=[0;1];
conditional.central.rho=[1 -1i;1i 1]/2;
for method={'cce','mecce','noise'}
    conditional.method=method{1};
    if strcmp(method{1},'mecce'), conditional.jumps=struct([]); end
    if strcmp(method{1},'noise')
        conditional.quantity='autocorrelation';
        conditional.noise=struct('terms',struct('support',2,'operator',diag([1 -1])/2),...
                                 'connected',true,'ordering','ordered','branch','alpha');
    end
    rejected=false;
    try
        cce(transverse,conditional);
    catch exception
        rejected=strcmp(exception.identifier,'cce:pure_dephasing');
    end
    result=test_true(result,['Conditional transition rejection: ' method{1}],rejected,...
                     'Discarding a nonzero central transition requires generalized dynamics.');

    % Exercise the same invariant-subspace guard through interlaced recursion
    if ~strcmp(method{1},'noise')
        ip=conditional; ip.bath_state=struct('mode','interlaced','basis','uniform',...
                                            'probabilities','uniform','samples',1,'seed',19);
        ip.interlaced=struct('condition_sites',{{[2 3]}},'reference_state',1,'max_states',4);
        ip.averaging='after_expansion'; rejected=false;
        try
            cce(transverse,ip);
        catch exception
            rejected=strcmp(exception.identifier,'cce:pure_dephasing');
        end
        result=test_true(result,['Interlaced transition rejection: ' method{1}],rejected,...
                         'Conditioning must not bypass the actual-Hamiltonian invariance check.');
    end
    if isfield(conditional,'jumps'), conditional=rmfield(conditional,'jumps'); end
end

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
    difference=Inf;
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

% Reject the obsolete control container through preparation and every evaluator
q.method='cce'; q.quantity='coherence'; q.sequence=[];
q=rmfield(q,{'noise','deer'});
q.controls=struct('bandwidth',1); canonical=cce_prepare(ss,rmfield(q,'controls'));
clean=rmfield(q,'controls'); ensemble=cce_states(canonical,clean); context=ensemble.contexts(1);
cluster=struct('sites',2,'physical_ids',2); schedule=cce_sequence(canonical,clean);
calls={@()cce_prepare(ss,q),@()cce(ss,q),...
       @()cce_cluster(canonical,cluster,context,schedule,q),...
       @()cce_master(canonical,cluster,context,schedule,q),...
       @()cce_noise(canonical,cluster,context,schedule,q),...
       @()cce_interlace(canonical,cluster,context,schedule,q),...
       @()cce_controls(canonical,[],[],schedule,q)};
for n=1:numel(calls)
    rejected=false;
    try
        calls{n}();
    catch exception
        rejected=contains(exception.message,'parameters.controls');
    end
    result=test_true(result,sprintf('Unsupported controls rejected path %d',n),rejected,...
                     'An accepted but unused control container must never change physics silently.');
end

% Audit ignored nested options and inactive driver branches
base=rmfield(q,{'controls','jumps'}); base.central=rmfield(base.central,'rho');
cases={base,base,base,base,base,base,base,base,base,base,base};
cases{1}.numerics.propagator='krylov'; cases{2}.numerics.gpu=true;
cases{3}.bath_state.temperature=300; cases{4}.bath_state.samples=2;
cases{5}.graph.cutoff=1; cases{6}.sampling=struct();
cases{7}.normalized=true; cases{8}.spectrum=struct('omega',0);
cases{9}.jumps=struct('support',2,'operator',eye(2),'rate',1);
cases{10}.sequence=struct('name','hahn','bandwidth',1);
cases{11}.sequence=struct('type','ideal','target',1,'placement','fraction',...
                          'at',0.5,'axis','x','angle',pi,'duration',1);
cases{12}=base; cases{12}.numerics.reduction='direct';
cases{13}=base; cases{13}.method='gcce'; cases{13}.central.rho=eye(2)/2;
cases{13}.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',0.1);
cases{14}=base; cases{14}.analysis=struct('fit',false,'model','exponential');
cases{15}=base; cases{15}.method='noise'; cases{15}.quantity='autocorrelation';
cases{15}.noise=struct('connected',false,'ordering','ordered',...
                      'terms',struct('support',2,'operator',eye(2),'weight',2));
cases{16}=cases{15}; cases{16}.noise.terms=rmfield(cases{16}.noise.terms,'weight');
cases{16}.jumps=struct('support',2,'operator',eye(2),'rate',1);
cases{17}=base; cases{17}.method='mecce';
cases{17}.jumps=struct('support',2,'operator',eye(2),'rate',1,'type','local');
labels={'propagator','gpu','temperature','samples','cutoff','sampling',...
        'normalized','spectrum','jumps','bandwidth','duration','reduction','second_order','fit','weight','jumps','type'};
for n=1:numel(cases)
    rejected=false;
    try
        cce(ss,cases{n});
    catch exception
        rejected=contains(exception.message,labels{n});
    end
    result=test_true(result,['Ignored option rejected: ' labels{n}],rejected,...
                     'Unsupported nested options must not silently leave a successful result unchanged.');
end

% Honor explicit raw suppression without changing the computed observable
retained=cce(ss,base); base.retain_raw=false; suppressed=cce(ss,base);
result=test_true(result,'Explicit raw suppression',isempty(suppressed.raw),...
                 'An accepted retain_raw=false request must not retain the raw cluster arrays.');
result=test_close(result,'Raw retention preserves physics',suppressed.value,retained.value,...
                  0,0,'Suppressing retained diagnostics must not change the observable.');

% Honor optional postprocessing through the public driver rather than dropping it
post=base; post=rmfield(post,'retain_raw'); post.method='noise'; post.quantity='autocorrelation';
post.noise=struct('connected',false,'ordering','ordered',...
                  'terms',struct('support',2,'operator',eye(2)));
post.spectrum=struct('frequency_grid',[0 1],'convention','exp(+i*omega*t)',...
                     'window','none','weights','trapezoid','sidedness','two-sided',...
                     'time_domain','positive','stationary',true);
post.analysis=struct('fit',false); actual=cce(ss,post);
reference=cce_spectrum(actual.times,actual.value,post);
result=test_true(result,'Driver retains requested postprocessing',...
                 isfield(actual,'analysis')&&isfield(actual,'spectrum'),...
                 'Accepted analysis and spectrum requests must produce their documented helper outputs.');
if isfield(actual,'spectrum'), difference=max(abs(actual.spectrum.value-reference.value),[],'all');
else, difference=Inf; end
result=test_close(result,'Driver spectrum helper parity',difference,0,...
                  0,0,'The driver must invoke the same explicit Fourier contract as cce_spectrum.');

% Retain physical mean fields while stripping their options during event recompilation
r=base; r=rmfield(r,intersect(fieldnames(r),{'noise','spectrum','deer','normalized'}));
r.quantity='coherence'; r.central=struct('alpha',[1;0],'beta',[0;1],'rho',ones(2)/2);
r.bath_state=struct('mode','mixed_product'); r.meanfield='static'; r.meanfield_radius=Inf;
r.times=[0 0.01]; r.order=2;
for method={'cce','mecce','megcce'}
    r.method=method{1};
    if ~strcmp(r.method,'cce'), r.jumps=struct([]); end
    for control=1:3

        % Cover named central pulses, DEER selection and explicit bath pulses
        r.sequence=struct('name','hahn','count',1,'axis','x','target',1);
        if control==2
            r.sequence=struct('type','ideal','target','deer','placement','fraction',...
                              'at',0.5,'axis','x','angle',pi);
            r.deer=struct('enabled',true,'detuning_center',0,'halfwidth',1,'sites',2);
        elseif control==3
            r.sequence=struct('type','ideal','target',2,'placement','fraction',...
                              'at',0.5,'axis','x','angle',pi);
        end
        if control~=2&&isfield(r,'deer'), r=rmfield(r,'deer'); end

        % Compare the explicit unbounded radius against the default unbounded field
        unlimited=rmfield(r,'meanfield_radius'); reference=cce_sequence(prepared,unlimited);
        passed=true; difference=Inf;
        try
            with_radius=cce_sequence(prepared,r); difference=double(~isequaln(with_radius,reference));
        catch exception
            if ~strcmp(exception.identifier,'cce:capability'), rethrow(exception); end
            passed=false;
        end
        name=sprintf('Static radius recursive controls %s %d',r.method,control);
        result=test_true(result,[name ' accepted'],passed,...
                         'Internal compilation must remove the radius when disabling mean fields.');
        result=test_close(result,[name ' parity'],difference,0,1e-12,1e-12,...
                          'An explicit infinite radius must equal the unbounded static-field default.');
    end
end

% Require actual geometry for every finite static-field radius path
result=radius_checks(result,ss,base);

fprintf('DEER exact-boundary cases passed: selected sites [3 4], [3 6], and empty.\n');

end

% Exercise all driver, direct-evaluator and event-recompilation paths
function result=radius_checks(result,spin_system,parameters)

% Use a polarized exterior bath so dropping its fields is not innocuous
parameters=rmfield(parameters,intersect(fieldnames(parameters),{'deer','noise'}));
parameters.central=struct('alpha',[1;0],'beta',[0;1],'rho',ones(2)/2);
parameters.bath_state=struct('mode','explicit_product',...
                             'rho',{{diag([1 0]),diag([1 0])}});
parameters.meanfield='static'; parameters.meanfield_radius=2;
parameters.order=1; parameters.sequence=[];
cluster=struct('id',1,'sites',2);
for method={'cce','gcce','mecce','megcce'}

    % Choose each implemented evaluator without changing the missing geometry
    parameters.method=method{1}; parameters.quantity='coherence';
    parameters.jumps=struct([]);
    model=cce_prepare(spin_system,parameters); ensemble=cce_states(model,parameters);
    context=ensemble.contexts(1); schedule=cce_sequence(model,parameters);

    % Distinguish missing central and exterior sites from irrelevant included sites
    for missing=[model.central_ids 3]
        partial=model; partial.coordinates_A=[0 0 0;1 0 0;3 0 0];
        partial.coordinates_A(missing,:)=NaN; rejected=false;
        try
            cce_cluster(partial,cluster,context,schedule,parameters);
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates')&&...
                     contains(err.message,'meanfield_radius');
        end
        result=test_true(result,sprintf('Finite radius missing site %d %s',missing,parameters.method),...
                         rejected,'Each central and exterior site needs finite coordinates.');
    end
    spatial=model; spatial.coordinates_A=[0 0 0;NaN NaN NaN;3 0 0];
    accepted=true;
    try
        if ismember(parameters.method,{'mecce','megcce'})
            cce_master(spatial,cluster,context,schedule,parameters);
        else
            cce_cluster(spatial,cluster,context,schedule,parameters);
        end
    catch
        accepted=false;
    end
    result=test_true(result,['Finite radius ignores included-site geometry ' parameters.method],...
                     accepted,'The radius does not inspect sites already included in the cluster.');

    % Driver and direct helper calls must report the same prerequisite
    master=ismember(parameters.method,{'mecce','megcce'});
    calls=cell(1,5+master); names=cell(size(calls));
    calls(1:3)={@()cce(spin_system,parameters),...
                @()cce_cluster(model,cluster,context,schedule,parameters),...
                @()cce_controls(model,cluster.sites,context,schedule,parameters)};
    names(1:3)={'driver','cluster','controls'};
    if master
        calls{6}=@()cce_master(model,cluster,context,schedule,parameters);
        names{6}='master';
    end
    % A controlled schedule must not bypass the pre-dispatch prerequisite
    driven=schedule; driven.controls=true;
    calls{4}=@()cce_cluster(model,cluster,context,driven,parameters);
    names{4}='controlled dispatch';

    % Conditional enumeration must use the same checked evaluators
    ip=parameters; ip.bath_state=struct('mode','interlaced','basis','uniform',...
                                      'probabilities','uniform','samples',1,'seed',7);
    ip.interlaced=struct('condition_sites',{{[2 3]}},'reference_state',1,'max_states',4);
    sampled=cce_states(model,ip);
    calls{5}=@()cce_interlace(model,cluster,sampled,schedule,ip);
    names{5}='interlace';
    for k=1:numel(calls)
        rejected=false;
        try
            calls{k}();
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates')&&...
                     contains(err.message,'meanfield_radius');
        end
        result=test_true(result,['Finite radius geometry ' parameters.method ' ' names{k}],...
                         rejected,'A finite mean-field radius requires central and exterior coordinates.');
    end
end

% Compile before a context exists, then validate geometry during propagation
parameters.method='gcce'; parameters.quantity='density'; parameters.times=[0 0.01];
parameters.sequence=struct('name','hahn','count',1,'axis','x','target',1);
parameters.meanfield_radius=2;
model=cce_prepare(spin_system,parameters); models={model,model};
models{2}.coordinates_A=[0 0 0;1 0 0;3 0 0];
for geometry=1:2
    compiled=false;
    try
        events=cce_controls(models{geometry},[],[],[],parameters);
        compiled=isscalar(events.schedule.events);
    catch
    end
    result=test_true(result,sprintf('Finite-radius pre-context compilation %d',geometry),...
                     compiled,'Control compilation requires no state or geometry; propagation does.');
end

% Exercise every route into control compilation with complete physical geometry
spatial=spin_system; spatial.inter.coordinates={[0 0 0],[1 0 0],[3 0 0]};
sequences={parameters.sequence,[],...
           struct('type','finite','target',1,'placement','absolute','at',0,...
                  'duration',0.005,'axis','x','amplitude',1,'detuning',0),...
           struct('type','ideal','target',2,'placement','fraction','at',0.5,...
                  'axis','z','angle',pi)};
names={'named','deer','finite','bath'};
for method={'cce','gcce','mecce','megcce'}
    for control=1:numel(sequences)

        % Finite drives and exterior bath rotations use the generalized evaluator
        if control>=3&&~strcmp(method{1},'gcce'), continue; end
        p=parameters; p.method=method{1}; p.quantity='coherence';
        p.sequence=sequences{control};
        if ismember(p.method,{'gcce','megcce'}), p.normalized=true; end
        if ismember(p.method,{'mecce','megcce'})
            p.jumps=struct('support',{},'operator',{},'rate',{});
        end
        if control==2
            p.deer=struct('enabled',true,'detuning_center',0,...
                          'halfwidth',1,'sites',2);
            p.sequence=sequences{4}; p.sequence.target='deer';
            if ~strcmp(p.method,'gcce')
                p.deer.halfwidth=0; p.deer.allow_empty=true;
            end
        end
        label=[p.method ' ' names{control}]; accepted=false; difference=Inf;
        try
            observed=cce(spatial,p); accepted=all(observed.valid,'all');
            difference=max(abs(observed.value(:)-1));
        catch err
            fprintf('Finite-radius %s: %s %s\n',label,err.identifier,err.message);
        end
        result=test_true(result,['Finite-radius control dispatch ' label],accepted,...
                         'Valid geometry and a finite static-field cutoff must survive compilation.');
        result=test_close(result,['Finite-radius identity control ' label],difference,0,...
                          1e-12,0,'An uncoupled zero-field register prepared along x retains unit coherence.');

        % The same controls must still reject unavailable propagation geometry
        rejected=false;
        try
            cce(spin_system,p);
        catch err
            rejected=strcmp(err.identifier,'cce:coordinates');
        end
        result=test_true(result,['Finite-radius propagation geometry ' label],rejected,...
                         'Deferring the context-dependent check must not remove the physical prerequisite.');
    end
end

% Keep invalid scalar radii from silently disabling or emptying the field
for radius={NaN,-1,-Inf,[1 2],1i,'2'}
    parameters.meanfield_radius=radius{1}; rejected=false;
    try
        cce_validate(parameters);
    catch err
        rejected=strcmp(err.identifier,'cce:parameters')&&contains(err.message,'meanfield_radius');
    end
    result=test_true(result,'Invalid finite-radius declaration',rejected,...
                     'meanfield_radius must be a nonnegative real scalar or positive Inf.');

    % Retain scalar validation even when compilation has no context
    rejected=false;
    try
        cce_controls(model,[],[],[],parameters);
    catch err
        rejected=strcmp(err.identifier,'cce:parameters')&&contains(err.message,'meanfield_radius');
    end
    result=test_true(result,'Invalid pre-context radius declaration',rejected,...
                     'Deferring geometry must not accept an invalid cutoff during compilation.');
end

end


