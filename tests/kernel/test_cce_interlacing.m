% Test: kernel/cce_interlacing | CCE interlacing
% Exact W12 conditioning-set and probability regression tests. Syntax:
%
%                    result=test_cce_interlacing()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_interlacing()

% Construct overlapping superclusters around a single bath spin
result=new_test_result('kernel/cce_interlacing','CCE projected interlacing',...
                       'Check unique configurations, conditional weights and shared denominator K.');

% Compile the overlapping-pair model and its conditional bath ensemble
sys.isotopes={'1H','1H','1H','1H'};
sys.magnet=0;
sys.disable={'hygiene'};
sys.parallel={'processes',1};
sys.output='hush';
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([0 0 4])/(2*pi);
inter.coupling.matrix{1,4}=diag([0 0 7])/(2*pi);
spin_system=create(sys,inter);
p.observed=1;
p.cce.bath_spins=2:4;
p.cce.approximation='generalized';
p.normalized=false;
p.cce.expansion='standard';
p.output='signal';
p.times=[0 0.2 0.4];
p.cce.order=2;
p.rho0=[1 1;1 1]/2;
p.transition.spins=1;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','complete','score','coupling_norm');
p.cce.bath_state=struct('mode','interlaced','samples',1,'seed',31,'basis','uniform',...
                        'probabilities',{{[0.7 0.3],[0.2 0.8],[0.6 0.4]}});
p.cce.meanfield='static';
p.cce.interlaced=struct('max_states',16,'reference_state',1);
p.cce.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
[~,details]=clexpn(spin_system,@free_evolution,p,'labframe');
p=details.parameters;
model=cce_compile(spin_system,p,'labframe');
model.clusters=details.clusters;
ensemble=details.ensemble;
schedule=struct('pulse_sequence',@free_evolution,'contract',details.contract);

pulse_sequence=@free_evolution;
evaluator=schedule.contract.evaluator;

% Refuse count-scaled storage before constructing conditioning assignments
q=p;
q.times=0;
q.cce.approximation='projected';
q.normalized=true;
q.cce.numerics.max_memory_bytes=1000;
q.cce.interlaced=struct('max_states',8,'reference_state',1,...
                        'condition_sites',{{[2 3 4]}},'raw_only',true);
probe=struct('id',uint64(0),'sites',[],'groups',[]);
refused=false;
try
    cce_interlace(model,probe,ensemble,schedule,q,evaluator);
catch exception
    refused=strcmp(exception.identifier,'cce_interlace:memory');
end
result=test_true(result,'Conditioning union memory preflight',refused,...
                 'The assignment union and sorting workspace must fit the byte budget, not just max_states.');

% Sampling is bounded before allocating configuration records or assignments
q=p;
q.cce.order=0;
q.cce.bath_state.samples=1024;
q.cce.numerics.max_memory_bytes=4096;
refused=false;
try
    clexpn(model.spin_system,@free_evolution,q,'labframe');
catch exception
    refused=strcmp(exception.identifier,'clexpn:memory');
end
result=test_true(result,'Sampled configuration memory preflight',refused,...
                 'Direct state compilation must reject samples whose storage exceeds the byte budget.');

% Direct propagation has the same output-storage contract as the scheduler
q=p;
q.times=linspace(0,1,1024);
q.cce.numerics.max_memory_bytes=4096;
refused=false;
try
    evaluator(model,probe,ensemble.spin_config(1),schedule,q);
catch exception
    refused=strcmp(exception.identifier,'clexpn:memory');
end
result=test_true(result,'Direct cluster output memory preflight',refused,...
                 'A short Hilbert space does not bound the number of retained time points.');
ensemble.spin_config(1).assignment=[1 1 1];
ensemble.spin_config(1).bath_rho={diag([1 0]),diag([1 0]),diag([1 0])};
p.cce.interlaced=struct('reference_state',1,'max_states',16,'condition_sites',{{[2 3],[2 4]}});
cluster=struct('id',uint64(1),'sites',2);
raw=cce_interlace(model,cluster,ensemble,schedule,p,evaluator);

% Reject reference states inconsistent with their declared product assignment
rejected=false;
supplied=ensemble;
supplied.spin_config(1).bath_rho{1}=eye(2)/2;
try
    cce_interlace(model,cluster,supplied,schedule,p,evaluator);
catch exception
    rejected=strcmp(exception.identifier,'cce:spin_config');
end
result=test_true(result,'Interlaced reference density consistency',rejected,...
                 'Interlacing cannot silently replace a supplied density inconsistent with its assignment.');

% Reject a conflicting spin_config assignment rather than replacing it silently
supplied=ensemble;
supplied.spin_config(1).assignment=[2 1 1];
rejected=false;
try
    cce_interlace(model,cluster,supplied,schedule,p,evaluator);
catch exception
    rejected=strcmp(exception.identifier,'cce:spin_config');
end
result=test_true(result,'Interlaced reference assignment consistency',rejected,...
                 'A supplied spin_config assignment must agree with the ensemble reference.');

% Reject an unsupported reference projection before averaging it away
supplied=ensemble;
supplied.spin_config(1).projection=eye(2);
rejected=false;
try
    cce_interlace(model,cluster,supplied,schedule,p,evaluator);
catch exception
    rejected=strcmp(exception.identifier,'cce_subsystem:configuration');
end
result=test_true(result,'Interlaced rejects spin_config projection',rejected,...
                 'Interlacing must reject unsupported supplied projection fields.');

% Propagate experiment-owned observed preparation through the conditional bath
q=p;
q.cce.approximation='generalized';
q.normalized=false;
q.output='density';
q.times=0;
supplied=ensemble;
q.rho0=diag([0.2 0.8]);
initial=cce_interlace(model,cluster,supplied,schedule,q,evaluator);
expected=q.rho0;
result=test_close(result,'Interlaced supplied observed state',initial.value,expected(:).',...
                  1e-12,1e-12,'The experiment-owned observed density survives the projected bath average.');

% Keep the empty-cluster projected baseline valid for an empty family
q=p;
q.cce.interlaced=rmfield(q.cce.interlaced,'condition_sites');
q.cce.order=0;
empty_cluster=struct('id',uint64(0),'sites',[]);
baseline_ok=true;
baseline_error=Inf;
try
    empty_model=model;
    empty_model.clusters=model.clusters([]);
    baseline=cce_interlace(empty_model,empty_cluster,ensemble,schedule,q,evaluator);
    reference=ensemble.spin_config(1);
    reference.bath_rho={diag([1 0]),diag([1 0]),diag([1 0])};
    direct=evaluator(model,empty_cluster,reference,schedule,q);
    baseline_error=norm(baseline.value-direct.value);
catch exception
    baseline_ok=false;
    fprintf('empty interlaced baseline exception %s\n',exception.identifier);
end
result=test_true(result,'Empty-family interlaced baseline',baseline_ok,...
                 'An empty family still has a well-defined baseline at the reference assignment.');
result=test_close(result,'Empty-family baseline parity',baseline_error,0,1e-12,1e-12,...
                  'No conditioning sites means exactly the unconditioned reference baseline.');
for method={'cce','gcce','mecce','megcce'}
    q.cce.order=0;
    q.cce.approximation='projected';
    q.normalized=true;
    if ismember(method{1},{'gcce','megcce'})
        q.cce.approximation='generalized';
        q.normalized=false;
    end
    q.cce.basis.formalism='zeeman-hilb';
    if ismember(method{1},{'mecce','megcce'})
        q.cce.basis.formalism='zeeman-liouv';
    end
    q.bath_jumps=struct([]);
    if strcmp(q.cce.approximation,'generalized')
        q.rho0=[1 1;1 1]/2;
    elseif isfield(q,'rho0')
        q=rmfield(q,'rho0');
    end
    try
        [driver,~]=clexpn(model.spin_system,pulse_sequence,q,'labframe');
        driver_ok=all(isfinite(driver(:)));
    catch exception
        driver_ok=false;
        fprintf('order-zero %s driver exception %s\n',method{1},exception.identifier);
    end
    result=test_true(result,['Order-zero interlaced ' method{1} ' driver'],driver_ok,...
                     'Every interlaced evaluator must retain its baseline with no clusters.');
end
fprintf('empty interlaced baseline error %.3g\n',baseline_error);

% Restore the overlapping-supercluster reference assignments
expected=[1 1 1;1 1 2;1 2 1;2 1 1;2 1 2;2 2 1];
weights=[0.7*0.2*0.6;0.7*0.2*0.4;0.7*0.8*0.6;...
         0.3*0.2*0.6;0.3*0.2*0.4;0.3*0.8*0.6];
weights=weights/sum(weights);
result=test_true(result,'Overlapping cylinders yield six unique configurations, not eight occurrences',isequal(raw.assignments,expected),'Overlapping cylinders yield six unique configurations, not eight occurrences');
result=test_close(result,'Projected probabilities are physical product weights',raw.weights,weights,1e-15,0,'Projected probabilities are physical product weights');

% Independently enumerate the full bath and accept the union predicate
value=zeros(size(raw.value));
mass=0;
duplicated=zeros(size(raw.value));
occurrence_mass=0;
for a=1:2
    for b=1:2
        for c=1:2
            multiplicity=(c==1)+(b==1);
            if multiplicity==0
                continue;
            end
            spin_config=ensemble.spin_config(1);
            assignment=[a b c];
            weight=1;
            for j=1:3
                rho=zeros(2);
                rho(assignment(j),assignment(j))=1;
                spin_config.bath_rho{j}=rho;
                weight=weight*p.cce.bath_state.probabilities{j}(assignment(j));
            end
            item=evaluator(model,cluster,spin_config,schedule,p);
            value=value+weight*item.value;
            mass=mass+weight;
            duplicated=duplicated+multiplicity*weight*item.value;
            occurrence_mass=occurrence_mass+multiplicity*weight;
        end
    end
end
value=value/mass;
duplicated=duplicated/occurrence_mass;
difference=norm(value-duplicated);
result=test_close(result,'Projected raw average agrees with independent full-state enumeration',raw.value,value,1e-12,0,'Projected raw average agrees with independent full-state enumeration');
result=test_true(result,'Duplicate-occurrence weighting gives a detectably different answer',difference>1e-4,'Duplicate-occurrence weighting gives a detectably different answer');

% Preserve the largest numerical residual across every projected assignment
probe=p;
probe.cce.interlaced.raw_only=true;
probe.cce.interlaced.condition_sites={[2 3 4]};
probe.cce.bath_state.basis={eye(2),[0 1;1 0],eye(2)};
probe_model=model;
probe_model.zeeman{2}(3,:)=[0.731 0.473 0.122];
sample=probe;
sample.cce.interlaced.raw_only=false;
[~,probe_details]=clexpn(probe_model.spin_system,@free_evolution,sample,'labframe');
probe_ensemble=probe_details.ensemble;
measured=cce_interlace(probe_model,cluster,probe_ensemble,schedule,probe,evaluator);
residual=0;
orthogonality=0;
for k=1:size(measured.assignments,1)
    spin_config=probe_ensemble.spin_config(1);
    for j=1:3
        vector=probe.cce.bath_state.basis{j}(:,measured.assignments(k,j));
        spin_config.bath_rho{j}=vector*vector';
    end
    item=evaluator(probe_model,cluster,spin_config,schedule,probe);
    residual=max(residual,max(item.eigen_residual));
    orthogonality=max(orthogonality,max(item.orthogonality_residual));
end
result=test_close(result,'Interlaced worst residual retained',max(measured.eigen_residual),residual,...
                  0,0,'Ensemble diagnostics must not report only the first assignment residual.');
result=test_close(result,'Interlaced worst orthogonality retained',max(measured.orthogonality_residual),...
                  orthogonality,0,0,'Residual metadata follows the worst constituent solve, not its weight.');

% Denominator subclusters are conditioned on exactly the numerator distribution
empty=struct('id',uint64(0),'sites',[]);
baseline=cce_interlace(model,empty,ensemble,schedule,p,evaluator);
result=test_true(result,...
                 'Explicit shared K is identical in numerator and empty denominator',isequal(raw.assignments,baseline.assignments)&&isequal(raw.weights,baseline.weights),...
                 'Explicit shared K is identical in numerator and empty denominator');
p.cce.interlaced.max_states=7;
refused=false;
try
    cce_interlace(model,cluster,ensemble,schedule,p,evaluator);
catch exception
    refused=strcmp(exception.identifier,'cce_interlace:budget');
end
result=test_true(result,'Enumeration work exceeding the declared budget is refused',refused,'Enumeration work exceeding the declared budget is refused');
p.cce.interlaced.max_states=16;
p.cce.interlaced=rmfield(p.cce.interlaced,'condition_sites');
automatic=cce_interlace(model,cluster,ensemble,schedule,p,evaluator);
result=test_true(result,'Included-supercluster family generates the same K automatically',isequal(automatic.assignments,expected),'Included-supercluster family generates the same K automatically');

% At full order every conditioning cylinder is the full bath and W12 is exact
p.cce.order=3;
p.cce.interlaced.max_states=32;
[integrated,~]=clexpn(model.spin_system,pulse_sequence,p,'labframe');
mixture=p;
mixture.cce.bath_state=struct('mode','explicit_product',...
                              'rho',{{diag([0.7 0.3]),diag([0.2 0.8]),diag([0.6 0.4])}});
mixture.cce.meanfield='none';
mixture.cce=rmfield(mixture.cce,'interlaced');
[exact,~]=clexpn(model.spin_system,pulse_sequence,mixture,'labframe');
result=test_close(result,'Full-order W12 equals exact mixed bath',integrated,exact,...
                  1e-11,0,'Shared-K recursion telescopes at full bath order');

% Full-bath interlacing bypasses singular lower-cluster ratios
sys.isotopes={'13C','13C','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 2]);
inter.coupling.matrix{1,3}=zeros(3);
inter.coupling.matrix{2,3}=diag([8 0 0]);
spin_system=create(sys,inter);
spin_system.inter.coupling.matrix{1,2}=diag([0 0 2]);
spin_system.inter.coupling.matrix{1,3}=zeros(3);
spin_system.inter.coupling.matrix{2,3}=diag([8 0 0]);
q=p;
q.cce.approximation='projected';
q.normalized=true;
q=rmfield(q,'rho0');
q.cce.bath_spins=[2 3];
q.cce.order=2;
q.times=[0 pi/2];
q.cce.bath_state.probabilities='uniform';
reference=q;
reference.cce.bath_state=struct('mode','mixed_product');
reference.cce=rmfield(reference.cce,'interlaced');
reference.cce.meanfield='none';
[exact,~]=clexpn(spin_system,pulse_sequence,reference,'labframe');
[full_bath,full_bath_details]=clexpn(spin_system,pulse_sequence,q,'labframe');

% Verify the raw full-bath value against independent projected exponentials
x=[0 1;1 0]/2;
z=diag([1 -1])/2;
ha=kron(z,eye(2))+8*kron(x,x);
hb=-kron(z,eye(2))+8*kron(x,x);
dense=trace(expm(1i*hb*pi/2)*expm(-1i*ha*pi/2))/4;
result=test_close(result,'Interlaced independent dense reference',full_bath(end),dense,...
                  1e-12,1e-12,'The full-bath bypass must reproduce a separately constructed four-dimensional trace.');
result=test_close(result,'Interlaced full-bath zero baseline',full_bath,exact,...
                  1e-12,1e-12,'The exact projected full-bath average survives a zero empty-cluster baseline.');
result=test_true(result,'Interlaced full-bath validity',all(full_bath_details.reduction.valid),...
                 'A singular auxiliary ratio must not invalidate an exact full-bath value.');

% Preserve the same bypass for raw generalised observed signal
q.cce.approximation='generalized';
q.normalized=false;
q.rho0=p.rho0;
[generalised,generalised_details]=clexpn(spin_system,pulse_sequence,q,'labframe');
result=test_close(result,'Generalized interlaced full-bath zero baseline',generalised,...
                  full_bath/2,1e-12,1e-12,...
                  'Raw generalized signal differs only by its initial amplitude of one half.');
fprintf('R1242 interlaced %.16g; exact %.16g.\n',real(full_bath(end)),real(exact(end)));

% Require the common result schema and diagnostics on the interlaced path
reference=q;
reference.cce.bath_state=struct('mode','mixed_product');
reference.cce=rmfield(reference.cce,'interlaced');
[standard,standard_details]=clexpn(spin_system,pulse_sequence,reference,'labframe');
result=test_true(result,'Interlaced common result fields',...
                 isnumeric(generalised)&&isnumeric(standard)&&isequal(size(generalised),size(standard)),...
                 'Every native context returns the same numeric time-by-channel contract.');
result=test_true(result,'Interlaced common diagnostic fields',...
                 isequal(sort(fieldnames(generalised_details)),sort(fieldnames(standard_details))),...
                 'Interlaced results must carry numerical coverage and budget diagnostics.');
q.cce.convergence=struct('ladders',struct('field','cce.order','values',{{2,2}}),...
                         'max_work',10,'atol',1e-8,'rtol',1e-8,'checks',{{}});
try
    [~,checked_details]=cce_converge(spin_system,q,pulse_sequence,'labframe');
    converged=strcmp(checked_details.convergence,'converged');
catch
    converged=false;
end
result=test_true(result,'Interlaced convergence diagnostics',converged,...
                 'A complete exact interlaced calculation must be usable by cce_converge.');
result=test_true(result,'Exact interlaced lower orders unavailable',...
                 all(isnan(generalised_details.reduction.by_order(:,:,1:end-1)),'all')&&...
                 all(~generalised_details.reduction.by_order_valid(:,:,1:end-1),'all'),...
                 'Bypassed lower orders must not masquerade as evaluated exact prefixes.');
result=test_close(result,'Exact interlaced final prefix',generalised_details.reduction.by_order(:,:,end),...
                  generalised,0,0,'The final prefix must contain the exact raw average.');

% Match ordinary master evaluation for recursive and full-bath interlacing
sys.isotopes={'1H','1H','1H'};
sys.magnet=0;
sys.output='hush';
inter.zeeman.scalar={0,0,0};
inter.coupling=struct('scalar',{cell(3)});
spin_system=create(sys,inter);
q=p;
q.cce.bath_spins=[2 3];
q.times=[0 0.3 0.8];
q.cce.meanfield='none';
q.cce.bath_state.probabilities='uniform';
q.cce.interlaced=struct('max_states',4);
q.bath_jumps=struct('support',1,'operator',diag([1 -1]),'rate',0.7);
for method={'mecce','megcce'}
    q.cce.approximation='projected';
    q.normalized=true;
    if strcmp(method{1},'megcce')
        q.cce.approximation='generalized';
        q.normalized=false;
        q.rho0=p.rho0;
    elseif isfield(q,'rho0')
        q=rmfield(q,'rho0');
    end
    q.cce.basis.formalism='zeeman-liouv';
    for order=1:2
        q.cce.order=order;
        reference=q;
        reference.cce.bath_state=struct('mode','mixed_product');
        reference.cce=rmfield(reference.cce,'interlaced');
        [expected,~]=clexpn(spin_system,pulse_sequence,reference,'labframe');
        matched=false;
        difference_me=Inf;
        try
            [actual,actual_details]=clexpn(spin_system,pulse_sequence,q,'labframe');
            difference_me=max(abs(actual-expected));
            matched=all(actual_details.reduction.valid)&&difference_me<2e-11;
        catch exception
            fprintf('Interlaced %s order %d: %s\n',method{1},order,exception.message);
        end
        result=test_true(result,sprintf('Interlaced %s order %d dispatch',method{1},order),matched,...
                         'Both recursive subset averages and full-bath shortcuts must retain Lindblad jumps.');
        scale=1;
        if strcmp(method{1},'megcce')
            scale=0.5;
        end
        result=test_close(result,sprintf('Analytic %s dephasing reference',method{1}),expected,...
                          scale*exp(-1.4*q.times(:)),2e-11,2e-11,...
                          'A observed longitudinal jump damps signal at twice its supplied rate.');
        fprintf('R5575 %s order %d error %.3g.\n',method{1},order,difference_me);
    end
end

fprintf('interlacing configurations %d; raw enumeration error %.3g; duplicate-weight difference %.6g; mass %.6g\n',...
        size(raw.assignments,1),norm(raw.value-value),difference,mass);

% Preserve driven raw schemas through both truncated and full-bath interlacing
sys.isotopes={'1H','1H','1H'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([1 3 4])/(2*pi);
spin_system=create(sys,inter);
q=p;
q.cce.approximation='generalized';
q.normalized=false;
q.cce.bath_spins=[2 3];
q.cce.meanfield='none';
q.cce.bath_state.probabilities='uniform';
q.times=[0 0.2 0.5];
q.cce.interlaced=struct('max_states',4);
driven_ok=true;
driven_error=0;
finite=struct('spins',1,'time',0.05,'duration',0.05,...
              'axis','x','amplitude',2,'detuning',0);
for sequence={@hahn_echo,@finite_drive}
    pulse_sequence=sequence{1};
    if isequal(pulse_sequence,@finite_drive)
        q.drive=finite;
        q.cce.numerics=rmfield(q.cce.numerics,'propagator');
    else
        q.cce.numerics.propagator='expm';
    end
    for order=1:2
        q.cce.order=order;
        reference=q;
        reference.cce.bath_state=struct('mode','mixed_product');
        reference.cce=rmfield(reference.cce,'interlaced');
        [expected,~]=clexpn(spin_system,pulse_sequence,reference,'labframe');
        try
            [actual,actual_details]=clexpn(spin_system,pulse_sequence,q,'labframe');
            driven_error=max(driven_error,norm(actual-expected,'fro'));
            driven_ok=driven_ok&&all(isfinite(actual),'all')&&...
                      all(isfinite(actual_details.reduction.abs_error_bound),'all');
        catch exception
            driven_ok=false;
            fprintf('Driven interlacing exception %s\n',exception.identifier);
        end
    end
end
result=test_true(result,'Driven interlaced raw schemas',driven_ok,...
                 'Named and finite controls must not require a static-only residual field.');
result=test_close(result,'Driven interlaced mixed-state equivalence',driven_error,0,2e-12,0,...
                  'Explicit uniform assignments and mixed local states give identical cluster averages.');
fprintf('R6 driven interlacing accepted %d; maximum error %.3g\n',driven_ok,driven_error);

end


