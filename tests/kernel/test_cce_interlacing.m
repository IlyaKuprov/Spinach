% Exact W12 conditioning-set and probability regression tests.
% Syntax: result=test_cce_interlacing()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_interlacing()

% Construct overlapping superclusters around a single bath spin
result=new_test_result('kernel/cce_interlacing','CCE conditional interlacing',...
                       'Check unique configurations, conditional weights and shared denominator K.');
[model,p]=fixture(); ensemble=cce_states(model,p);
ensemble.assignments(1,:)=[1 1 1]; schedule=cce_sequence(model,p);
p.interlaced=struct('reference_state',1,'max_states',16,'condition_sites',{{[2 3],[2 4]}});
cluster=struct('id',uint64(1),'sites',2); raw=cce_interlace(model,cluster,ensemble,schedule,p);
% Keep the empty-cluster conditional baseline valid for an empty family
q=p; q.interlaced=rmfield(q.interlaced,'condition_sites'); q.order=0;
empty_cluster=struct('id',uint64(0),'sites',[]); baseline_ok=true; baseline_error=Inf;
try
    baseline=cce_interlace(model,empty_cluster,ensemble,schedule,q);
    reference=ensemble.contexts(1); reference.assignment=ensemble.assignments(1,:);
    reference.bath_rho={diag([1 0]),diag([1 0]),diag([1 0])};
    reference.expectations=repmat([0 0 0.5],3,1);
    direct=cce_cluster(model,empty_cluster,reference,schedule,q);
    baseline_error=norm(baseline.value-direct.value);
catch exception
    baseline_ok=false; fprintf('empty interlaced baseline exception %s\n',exception.identifier);
end
result=test_true(result,'Empty-family interlaced baseline',baseline_ok,...
                 'An empty family still has a well-defined baseline at the reference assignment.');
result=test_close(result,'Empty-family baseline parity',baseline_error,0,1e-12,1e-12,...
                  'No conditioning sites means exactly the unconditioned reference baseline.');
for method={'cce','gcce','mecce','megcce'}
    q.order=0; q.method=method{1}; q.jumps=struct([]);
    q.central.rho=[1 1;1 1]/2;
    try
        driver=cce(model.spin_system,q); driver_ok=all(isfinite(driver.value(:)));
    catch exception
        driver_ok=false; fprintf('order-zero %s driver exception %s\n',method{1},exception.identifier);
    end
    result=test_true(result,['Order-zero interlaced ' method{1} ' driver'],driver_ok,...
                     'Every interlaced evaluator must retain its baseline with no clusters.');
end
fprintf('empty interlaced baseline error %.3g\n',baseline_error);

% Restore the overlapping-supercluster reference assignments
expected=[1 1 1;1 1 2;1 2 1;2 1 1;2 1 2;2 2 1];
weights=[0.7*0.2*0.6;0.7*0.2*0.4;0.7*0.8*0.6;...
         0.3*0.2*0.6;0.3*0.2*0.4;0.3*0.8*0.6]; weights=weights/sum(weights);
result=test_true(result,'Overlapping cylinders yield six unique configurations, not eight occurrences',isequal(raw.assignments,expected),'Overlapping cylinders yield six unique configurations, not eight occurrences');
result=test_close(result,'Conditional probabilities are physical product weights',raw.weights,weights,1e-15,0,'Conditional probabilities are physical product weights');

% Independently enumerate the full bath and accept the union predicate
value=zeros(size(raw.value)); mass=0; duplicated=zeros(size(raw.value)); occurrence_mass=0;
for a=1:2
    for b=1:2
        for c=1:2
            multiplicity=(c==1)+(b==1);
            if multiplicity==0, continue; end
            context=ensemble.contexts(1); assignment=[a b c]; weight=1;
            for j=1:3
                rho=zeros(2); rho(assignment(j),assignment(j))=1; context.bath_rho{j}=rho;
                context.expectations(j,:)=[0 0 1.5-assignment(j)];
                weight=weight*p.bath_state.probabilities{j}(assignment(j));
            end
            item=cce_cluster(model,cluster,context,schedule,p);
            value=value+weight*item.value; mass=mass+weight;
            duplicated=duplicated+multiplicity*weight*item.value; occurrence_mass=occurrence_mass+multiplicity*weight;
        end
    end
end
value=value/mass; duplicated=duplicated/occurrence_mass; difference=norm(value-duplicated);
result=test_close(result,'Conditional raw average agrees with independent full-state enumeration',raw.value,value,1e-12,0,'Conditional raw average agrees with independent full-state enumeration');
result=test_true(result,'Duplicate-occurrence weighting gives a detectably different answer',difference>1e-4,'Duplicate-occurrence weighting gives a detectably different answer');

% Preserve the largest numerical residual across every conditional assignment
probe=p; probe.interlaced.raw_only=true; probe.interlaced.condition_sites={[2 3 4]};
probe.bath_state.basis={eye(2),[0 1;1 0],eye(2)};
probe_model=model; probe_model.zeeman{2}(3,:)=[0.731 0.473 0.122];
probe_ensemble=cce_states(probe_model,probe);
measured=cce_interlace(probe_model,cluster,probe_ensemble,schedule,probe);
residual=0; orthogonality=0;
for k=1:size(measured.assignments,1)
    context=probe_ensemble.contexts(1);
    for j=1:3
        vector=probe.bath_state.basis{j}(:,measured.assignments(k,j));
        context.bath_rho{j}=vector*vector';
        context.expectations(j,:)=[0 0 trace(context.bath_rho{j}*diag([0.5 -0.5]))];
    end
    item=cce_cluster(probe_model,cluster,context,schedule,probe);
    residual=max(residual,max(item.eigen_residual));
    orthogonality=max(orthogonality,max(item.orthogonality_residual));
end
result=test_close(result,'Interlaced worst residual retained',max(measured.eigen_residual),residual,...
                  0,0,'Ensemble diagnostics must not report only the first assignment residual.');
result=test_close(result,'Interlaced worst orthogonality retained',max(measured.orthogonality_residual),...
                  orthogonality,0,0,'Residual metadata follows the worst constituent solve, not its weight.');

% Denominator subclusters are conditioned on exactly the numerator distribution
empty=struct('id',uint64(0),'sites',[]); baseline=cce_interlace(model,empty,ensemble,schedule,p);
result=test_true(result,...
                 'Explicit shared K is identical in numerator and empty denominator',isequal(raw.assignments,baseline.assignments)&&isequal(raw.weights,baseline.weights),...
                 'Explicit shared K is identical in numerator and empty denominator');
p.interlaced.max_states=7; refused=false;
try
    cce_interlace(model,cluster,ensemble,schedule,p); catch exception, refused=strcmp(exception.identifier,'cce_interlace:budget'); end
result=test_true(result,'Enumeration work exceeding the declared budget is refused',refused,'Enumeration work exceeding the declared budget is refused');
p.interlaced.max_states=16; p.interlaced=rmfield(p.interlaced,'condition_sites');
automatic=cce_interlace(model,cluster,ensemble,schedule,p);
result=test_true(result,'Included-supercluster family generates the same K automatically',isequal(automatic.assignments,expected),'Included-supercluster family generates the same K automatically');
% At full order every conditioning cylinder is the full bath and W12 is exact
p.order=3; p.interlaced.max_states=32; integrated=cce(model.spin_system,p);
mixture=p; mixture.bath_state=struct('mode','explicit_product',...
    'rho',{{diag([0.7 0.3]),diag([0.2 0.8]),diag([0.6 0.4])}});
mixture.meanfield='none'; mixture=rmfield(mixture,'interlaced');
exact=cce(model.spin_system,mixture);
result=test_close(result,'Full-order W12 equals exact mixed bath',integrated.value,exact.value,...
                  1e-11,0,'Shared-K recursion telescopes at full bath order');

% Full-bath interlacing bypasses singular lower-cluster ratios
sys.isotopes={'13C','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)}; inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 2]);
inter.coupling.matrix{1,3}=zeros(3);
inter.coupling.matrix{2,3}=diag([8 0 0]); spin_system=create(sys,inter);
spin_system.inter.coupling.matrix{1,2}=diag([0 0 2]);
spin_system.inter.coupling.matrix{1,3}=zeros(3);
spin_system.inter.coupling.matrix{2,3}=diag([8 0 0]);
q=p; q.method='cce'; q.bath_spins=[2 3]; q.order=2; q.times=[0 pi/2];
q.bath_state.probabilities='uniform';
reference=q; reference.bath_state=struct('mode','mixed_product'); reference=rmfield(reference,'interlaced'); reference.meanfield='none';
exact=cce(spin_system,reference); full_bath=cce(spin_system,q);

% Verify the raw full-bath value against independent conditional exponentials
x=[0 1;1 0]/2; z=diag([1 -1])/2;
ha=kron(z,eye(2))+8*kron(x,x); hb=-kron(z,eye(2))+8*kron(x,x);
dense=trace(expm(1i*hb*pi/2)*expm(-1i*ha*pi/2))/4;
result=test_close(result,'Interlaced independent dense reference',full_bath.value(end),dense,...
                  1e-12,1e-12,'The full-bath bypass must reproduce a separately constructed four-dimensional trace.');
result=test_close(result,'Interlaced full-bath zero baseline',full_bath.value,exact.value,...
                  1e-12,1e-12,'The exact conditional full-bath average survives a zero empty-cluster baseline.');
result=test_true(result,'Interlaced full-bath validity',all(full_bath.valid),...
                 'A singular auxiliary ratio must not invalidate an exact full-bath value.');
% Preserve the same bypass for raw generalized central coherence
q.method='gcce'; generalized=cce(spin_system,q);
result=test_close(result,'Generalized interlaced full-bath zero baseline',generalized.value,...
                  full_bath.value/2,1e-12,1e-12,...
                  'Raw generalized coherence differs only by its initial amplitude of one half.');
fprintf('R1242 interlaced %.16g; exact %.16g.\n',real(full_bath.value(end)),real(exact.value(end)));

% Require the common result schema and diagnostics on the interlaced path
reference=q; reference.bath_state=struct('mode','mixed_product'); reference=rmfield(reference,'interlaced'); standard=cce(spin_system,reference);
result=test_true(result,'Interlaced common result fields',...
                 isequal(sort(fieldnames(generalized)),sort(fieldnames(standard))),...
                 'All public driver branches must expose the same documented result fields.');
result=test_true(result,'Interlaced common diagnostic fields',...
                 isequal(sort(fieldnames(generalized.diagnostics)),sort(fieldnames(standard.diagnostics))),...
                 'Interlaced results must carry numerical coverage and budget diagnostics.');
q.convergence=struct('ladders',struct('field','order','values',{{2,2}}),...
                     'max_work',10,'atol',1e-8,'rtol',1e-8,'checks',{{}});
try
    checked=cce_converge(spin_system,q); converged=checked.diagnostics.acceptable;
catch
    converged=false;
end
result=test_true(result,'Interlaced convergence diagnostics',converged,...
                 'A complete exact interlaced calculation must be usable by cce_converge.');
result=test_true(result,'Exact interlaced lower orders unavailable',...
                 all(isnan(generalized.by_order(:,:,1:end-1)),'all')&&...
                 all(~generalized.by_order_valid(:,:,1:end-1),'all'),...
                 'Bypassed lower orders must not masquerade as evaluated exact prefixes.');
result=test_close(result,'Exact interlaced final prefix',generalized.by_order(:,:,end),...
                  generalized.value,0,0,'The final prefix must contain the exact raw average.');

% Match ordinary master evaluation for recursive and full-bath interlacing
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0,0,0}; inter.coupling=struct('scalar',{cell(3)});
spin_system=create(sys,inter);
q=p; q.bath_spins=[2 3]; q.times=[0 0.3 0.8]; q.meanfield='none';
q.bath_state.probabilities='uniform'; q.interlaced=struct('max_states',4);
q.jumps=struct('support',1,'operator',diag([1 -1]),'rate',0.7);
for method={'mecce','megcce'}
    q.method=method{1};
    for order=1:2
        q.order=order; reference=q; reference.bath_state=struct('mode','mixed_product'); reference=rmfield(reference,'interlaced');
        expected=cce(spin_system,reference); matched=false; difference_me=Inf;
        try
            actual=cce(spin_system,q); difference_me=max(abs(actual.value-expected.value));
            matched=all(actual.valid)&&difference_me<2e-11;
        catch exception
            fprintf('Interlaced %s order %d: %s\n',method{1},order,exception.message);
        end
        result=test_true(result,sprintf('Interlaced %s order %d dispatch',method{1},order),matched,...
                         'Both recursive subset averages and full-bath shortcuts must retain Lindblad jumps.');
        scale=1; if strcmp(method{1},'megcce'), scale=0.5; end
        result=test_close(result,sprintf('Analytic %s dephasing reference',method{1}),expected.value,...
                          scale*exp(-1.4*q.times(:)),2e-11,2e-11,...
                          'A central Pauli-z jump damps coherence at twice its supplied rate.');
        fprintf('R5575 %s order %d error %.3g.\n',method{1},order,difference_me);
    end
end

fprintf('interlacing configurations %d; raw enumeration error %.3g; duplicate-weight difference %.6g; mass %.6g\n',...
        size(raw.assignments,1),norm(raw.value-value),difference,mass);

% Preserve driven raw schemas through both truncated and full-bath interlacing
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'processes',1};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([1 3 4])/(2*pi); spin_system=create(sys,inter);
q=p; q.method='gcce'; q.bath_spins=[2 3]; q.meanfield='none';
q.bath_state.probabilities='uniform'; q.times=[0 0.2 0.5];
q.interlaced=struct('max_states',4); driven_ok=true; driven_error=0;
finite=struct('type','finite','target',1,'placement','absolute','at',0.05,...
              'duration',0.05,'axis','x','amplitude',2,'detuning',0);
for sequence={struct('name','hahn','count',1,'axis','x','target',1),finite}
    q.sequence=sequence{1};
    for order=1:2
        q.order=order; reference=q; reference.bath_state=struct('mode','mixed_product'); reference=rmfield(reference,'interlaced');
        expected=cce(spin_system,reference);
        try
            actual=cce(spin_system,q);
            driven_error=max(driven_error,norm(actual.value-expected.value,'fro'));
            driven_ok=driven_ok&&all(isfinite(actual.value),'all')&&...
                      all(isfinite(actual.abs_error_bound),'all');
        catch exception
            driven_ok=false; fprintf('Driven interlacing exception %s\n',exception.identifier);
        end
    end
end
result=test_true(result,'Driven interlaced raw schemas',driven_ok,...
                 'Named and finite controls must not require a static-only residual field.');
result=test_close(result,'Driven interlaced mixed-state equivalence',driven_error,0,2e-12,0,...
                  'Explicit uniform assignments and mixed local states give identical cluster averages.');
fprintf('R6 driven interlacing accepted %d; maximum error %.3g\n',driven_ok,driven_error);

end

% Build a nontrivial static-exterior model with two included pair superclusters
function [model,p]=fixture()
sys.isotopes={'1H','1H','1H','1H'}; sys.magnet=0; sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.output='hush';
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([0 0 4])/(2*pi);
inter.coupling.matrix{1,4}=diag([0 0 7])/(2*pi);
spin_system=create(sys,inter); p.central_spins=1; p.bath_spins=2:4;
p.method='gcce'; p.expansion='standard'; p.quantity='coherence'; p.times=[0 0.2 0.4]; p.order=2;
p.central.rho=[1 1;1 1]/2; p.central.alpha=[1;0]; p.central.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm');
p.bath_state=struct('mode','interlaced','samples',1,'seed',31,'basis','uniform',...
                    'probabilities',{{[0.7;0.3],[0.2;0.8],[0.6;0.4]}});
p.meanfield='static'; p.sequence=[];
p.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
model=cce_prepare(spin_system,p);
end


