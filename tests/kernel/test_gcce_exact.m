% Independent full-Hilbert gCCE validation with transverse exchange. Syntax:
%
%                         result=test_gcce_exact()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_exact()

% Compare every raw density with a separate dense Kronecker/expm calculation
result=new_test_result('kernel/gcce_exact','gCCE exact validation',...
                       'Four bath spins, complex pure/mixed states and physical Hahn pulses.');

% Construct deterministic observed-exchange and bath-bath interactions
nb=4;
sys.isotopes=repmat({'13C'},1,nb+1);
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,nb+1);
inter.coupling.matrix=cell(nb+1);
spin_system=create(sys,inter);
fields=[31 17 23 29 37];
tensors=cell(nb+1);
for j=1:nb+1
    spin_system.inter.zeeman.matrix{j}(3,3)=fields(j);
    for k=j+1:nb+1
        tensors{j,k}=[7 2 -3;-1 11 4;3 -2 13]/(j+k);
        spin_system.inter.coupling.matrix{j,k}=tensors{j,k}.';
    end
end

% Normalise independently scaled observables with scale-aware uncertainty
rho0=diag([1 0]);
q=struct('observed',1,'times',0,'output','observable',...
         'rho0',rho0,'normalized',true,...
         'coil',{{diag([1 0]),1e-20*diag([1 0]),1e20*diag([1 0])}});
q.cce=struct('bath_spins',[],'order',0,'approximation','generalized',...
             'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^29));
model=cce_compile(spin_system,q,'labframe');
observed_system=cce_subsystem(model,[],struct(),[]);
scaled=NaN(1,3);
accepted=true;
try
    [~,~,initial_signal,meta]=cce_observation(observed_system,q);
    scaled=initial_signal./meta.normalization_divisor;
catch exception
    accepted=false;
    fprintf('scaled observable exception %s\n',exception.identifier);
end
result=test_true(result,'Scaled observable normalization accepted',accepted,...
                 'Changing physical units cannot turn a well-resolved nonzero reference into zero.');
result=test_close(result,'Scaled observable normalization parity',scaled,ones(1,3),1e-13,1e-13,...
                  'All independently scaled projectors must have normalized initial value one.');
q.coil={1e20*diag([1 -1])};
q.rho0=diag([0.5 0.5]);
refused=false;
try
    cce_observation(observed_system,q);
catch exception
    refused=strcmp(exception.identifier,'cce_observation:reference');
end
result=test_true(result,'Scaled zero reference rejection',refused,...
                 'Scaling an operator must not hide a genuinely zero normalization reference.');
fprintf('scaled observable normalized values %s\n',mat2str(scaled,12));

p.observed=1;
p.cce.bath_spins=2:5;
p.cce.approximation='generalized';
p.cce.expansion='standard';
p.output='density';
p.cce.order=4;
v=[sqrt(0.3);1i*sqrt(0.7)];
p.rho0=v*v';
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','complete','score','coupling_norm');
p.cce.meanfield='none';
p.times=linspace(0,0.07,8);
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^29);
maximum=0;
trace_error=0;
hermiticity=0;
mineig=1;
for pure=[false true]
    states=repmat({eye(2)/2},1,nb);
    if pure
        for j=1:nb
            v=[cos(j/5);exp(1i*j/4)*sin(j/5)];
            states{j}=v*v';
        end
    end
    p.cce.bath_state=struct('mode','explicit_product','rho',{states});
    for hahn=[false true]

        pulse_sequence=@free_evolution;
        if hahn
            pulse_sequence=@hahn_echo;
            p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1i]);
        elseif isfield(p,'transition')
            p=rmfield(p,'transition');
        end
        [~,diagnostics]=clexpn(spin_system,pulse_sequence,p,'labframe');
        for record=diagnostics.raw{1}

            % Embed operators independently in descending physical spin order
            sites=[];
            if record.cluster_id>0
                sites=diagnostics.clusters(record.cluster_id).sites;
            end
            ids=[1 sites];
            n=numel(ids);
            dim=2^n;
            db=dim/2;
            local={[0 1;1 0]/2,[0 -1i;1i 0]/2,diag([1 -1])/2};
            ops=cell(n,3);
            H=zeros(dim);
            rho=p.rho0;
            for j=1:n
                for a=1:3
                    ops{j,a}=kron(eye(2^(j-1)),kron(local{a},eye(2^(n-j))));
                end
                H=H+fields(ids(j))*ops{j,3};
                if j>1
                    rho=kron(rho,states{ids(j)-1});
                end
            end
            for j=1:n
                for k=j+1:n
                    for a=1:3
                        for b=1:3
                            H=H+tensors{ids(j),ids(k)}(a,b)*ops{j,a}*ops{k,b};
                        end
                    end
                end
            end
            pulse=expm(-1i*pi*ops{1,2});
            reference=zeros(numel(p.times),4);
            for k=1:numel(p.times)
                t=p.times(k);
                u=expm(-1i*H*t);
                if hahn&&t>0
                    u=expm(-1i*H*t/2)*pulse*expm(-1i*H*t/2);
                end
                evolved=u*rho*u';
                observed=zeros(2);
                for a=1:2
                    for b=1:2
                        observed(a,b)=trace(evolved((a-1)*db+(1:db),(b-1)*db+(1:db)));
                    end
                end
                reference(k,:)=observed(:).';
            end
            discrepancy=max(abs(record.value-reference),[],'all');
            maximum=max(maximum,discrepancy);
            result=test_true(result,'raw full-observed-spin system exact',...
                             discrepancy<=max(record.abs_error_bound,[],'all')+100*dim^3*eps,...
                             'Independent analytic or dense reference.');

        end
        reduced=diagnostics.reduction;
        reduced.density=reshape(reduced.value.',2,2,[]);
        physical=cce_diagnostics(struct('observed_dimension',2),diagnostics.raw,reduced,...
                                 struct('times',p.times,'numerics',p.cce.numerics));
        trace_error=max(trace_error,max(physical.trace_error));
        hermiticity=max(hermiticity,max(physical.hermiticity_error));
        mineig=min(mineig,min(physical.minimum_eigenvalue));

    end
end
metrics.maximum_raw_error=maximum;
metrics.trace_error=trace_error;
metrics.hermiticity_error=hermiticity;
metrics.minimum_eigenvalue=mineig;

% Report the measured metrics
disp(metrics);

% Keep observable norm factors in output-column order on every evaluator path
sys.isotopes={'1H','1H','1H'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter=struct();
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([1 3 4])/(2*pi);
small=create(sys,inter);
q=p;
q.cce.bath_spins=[2 3];
q.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
q.output='observable';
q.rho0=diag([1 0]);
q.coil={eye(2),2*eye(2)};
q.cce.bath_state=struct('mode','mixed_product');
q.bath_jumps=struct('support',1,'operator',diag([1 -1]),'rate',0.7);
accepted=true;
value_error=0;
bound_error=0;
cases=0;
finite=struct('spins',1,'time',0.05,'duration',0.05,...
              'axis','x','amplitude',2,'detuning',0);
for method={'gcce','megcce'}
    q.cce.basis=struct('formalism','zeeman-hilb');
    if strcmp(method{1},'gcce')
        q.bath_jumps=struct([]);
    else
        q.bath_jumps=struct('support',1,'operator',diag([1 -1]),'rate',0.7);
        q.cce.basis.formalism='zeeman-liouv';
    end
    for sequence={@free_evolution,@hahn_echo,@finite_drive}
        pulse_sequence=sequence{1};
        if strcmp(method{1},'megcce')&&~isequal(pulse_sequence,@free_evolution)
            continue;
        end
        if isfield(q.cce.numerics,'propagator')
            q.cce.numerics=rmfield(q.cce.numerics,'propagator');
        end
        if isfield(q,'drive')
            q=rmfield(q,'drive');
        end
        if isequal(pulse_sequence,@hahn_echo)
            q.cce.numerics.propagator='expm';
        elseif isequal(pulse_sequence,@finite_drive)
            q.drive=finite;
        end
        for mode={'mixed_product','interlaced'}
            q.cce.bath_state=struct('mode',mode{1});
            if isfield(q.cce,'interlaced')
                q.cce=rmfield(q.cce,'interlaced');
            end
            if strcmp(mode{1},'interlaced')
                q.cce.bath_state.samples=1;
                q.cce.bath_state.seed=31;
                q.cce.bath_state.probabilities='uniform';
                q.cce.bath_state.basis='uniform';
                q.cce.interlaced.max_states=4;
            end
            for order=1:2
                q.cce.order=order;
                for nt=2:3
                    q.times=linspace(0,0.5,nt);
                    for normalized=[false true]
                        q.normalized=normalized;
                        q.coil={eye(2),2*eye(2)};
                        [expected,expected_bound]=column_evaluate(small,pulse_sequence,q);
                        q.coil=q.coil.';
                        cases=cases+1;
                        try
                            [actual,actual_bound]=column_evaluate(small,pulse_sequence,q);
                            value_error=max(value_error,norm(actual-expected,'fro'));
                            bound_error=max(bound_error,max(abs(actual_bound-...
                                            expected_bound)./expected_bound,[],'all'));
                        catch exception
                            accepted=false;
                            fprintf('Column observable exception %s\n',exception.identifier);
                        end
                    end
                end
            end
        end
    end
end
result=test_true(result,'Column observable cells accepted',accepted,...
                 'Column cells must work for static, driven, master-equation and interlaced calculations.');
result=test_close(result,'Column observable values',value_error,0,2e-13,0,...
                  'Cell orientation cannot change the order or values of observable columns.');
result=test_close(result,'Column observable uncertainty',bound_error,0,2e-13,0,...
                  'Norm factors belong to observable columns, not to time rows.');

% Mixed observable units must retain the same normalisation decisions in column cells
observe=q;
observe.coil={1e-20*eye(2),1e20*eye(2)};
observe.normalized=true;
[~,~,initial_signal,meta]=cce_observation(observed_system,observe);
row=initial_signal./meta.normalization_divisor;
observe.coil=observe.coil.';
column=NaN(size(row));
scaled_ok=true;
try
    [~,~,initial_signal,meta]=cce_observation(observed_system,observe);
    column=initial_signal./meta.normalization_divisor;
catch
    scaled_ok=false;
end
result=test_true(result,'Column scaled normalization accepted',scaled_ok,...
                 'Each reference must be compared only with its own operator scale.');
result=test_close(result,'Column scaled normalization values',column,row,1e-14,1e-14,...
                  'Independent unit changes commute with cell transposition.');

% An unresolved first observable must not be hidden by another valid column
observe.coil={diag([0 1]);eye(2)};
rejected=false;
try
    cce_observation(observed_system,observe);
catch exception
    rejected=strcmp(exception.identifier,'cce_observation:reference');
end
result=test_true(result,'Column zero normalization rejected',rejected,...
                 'Every observable reference must be resolved before normalization.');
fprintf('R6 column-cell cases %d; accepted %d; value error %.3g; relative bound error %.3g\n',...
        cases,accepted,value_error,bound_error);

% Unknown and representation-incompatible propagation preferences remain errors
for preference={'unknown','eig'}
    request=q;
    request.cce.numerics.propagator=preference{1};
    unsupported=false;
    try
        clexpn(small,@free_evolution,request,'labframe');
    catch exception
        unsupported=strcmp(exception.identifier,'clexpn:propagator');
    end
    result=test_true(result,'Native propagation preference refusal',unsupported,...
                     'A master-equation request must not silently discard its propagation preference.');
end

% A requested direct-segment method must not be ignored in another representation
request=observe;
request.cce.numerics.propagator='expm';
for projected=[true false]
    restricted=observed_system;
    generator=struct('left',0,'right',0);
    if ~projected
        restricted.formalism='zeeman-liouv';
        generator=zeros(4);
    end
    unsupported=false;
    try
        hahn_echo(restricted,request,generator,0,0);
    catch exception
        unsupported=strcmp(exception.identifier,'hahn_echo:propagator');
    end
    result=test_true(result,'Direct Hahn representation refusal',unsupported,...
                     'Direct segment exponentials require closed generalized Hilbert dynamics.');
end

end

% Keep raw interlaced averages distinct from reduced observable responses
function [value,bound]=column_evaluate(spin_system,pulse_sequence,parameters)

% The context owns configuration draws and the ordinary observable reduction
request=parameters;
if strcmp(parameters.cce.bath_state.mode,'interlaced')
    request.cce.order=0;
end
[value,diagnostics]=clexpn(spin_system,pulse_sequence,request,'labframe');
bound=diagnostics.reduction.abs_error_bound;
if ~strcmp(parameters.cce.bath_state.mode,'interlaced')
    return;
end

% Evaluate the same raw conditional average using the registered callback
[parameters,~]=cce_options(spin_system,parameters,pulse_sequence,'labframe');
model=cce_compile(spin_system,parameters,'labframe');
model.formalism=parameters.cce.basis.formalism;
model.clusters=diagnostics.clusters;
schedule=struct('pulse_sequence',pulse_sequence,'contract',diagnostics.contract);
parameters.cce.interlaced.raw_only=true;
parameters.cce.interlaced.reference_state=1;
cluster=struct('id',uint64(1),'sites',parameters.cce.bath_spins(1:parameters.cce.order));
raw=cce_interlace(model,cluster,diagnostics.ensemble,schedule,parameters,diagnostics.contract.evaluator);
value=raw.value;
bound=raw.abs_error_bound;

end


