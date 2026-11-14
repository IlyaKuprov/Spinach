% Independent full-Hilbert gCCE validation with transverse exchange. Syntax:
%
%                         result=test_gcce_exact()
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_exact()

% Compare every raw density with a separate dense Kronecker/expm calculation
result=new_test_result('kernel/gcce_exact','gCCE exact validation',...
                       'Four bath spins, complex pure/mixed states and physical Hahn pulses.');

% Normalize independently scaled observables with scale-aware uncertainty
rho0=diag([1 0]); q=struct('quantity','observable','central',struct('rho',rho0),...
                         'normalized',true,'operators',{{diag([1 0]),1e-20*diag([1 0]),1e20*diag([1 0])}});
scaled=NaN(1,3); accepted=true;
try
    scaled=cce_observe(rho0,q);
catch exception
    accepted=false; fprintf('scaled observable exception %s\n',exception.identifier);
end
result=test_true(result,'Scaled observable normalization accepted',accepted,...
                 'Changing physical units cannot turn a well-resolved nonzero reference into zero.');
result=test_close(result,'Scaled observable normalization parity',scaled,ones(1,3),1e-13,1e-13,...
                  'All independently scaled projectors must have normalized initial value one.');
q.operators={1e20*diag([1 -1])}; q.central.rho=diag([0.5 0.5]); refused=false;
try
    cce_observe(q.central.rho,q);
catch exception
    refused=strcmp(exception.identifier,'cce_observe:initial_zero');
end
result=test_true(result,'Scaled zero reference rejection',refused,...
                 'Scaling an operator must not hide a genuinely zero normalization reference.');
fprintf('scaled observable normalized values %s\n',mat2str(scaled,12));

% Construct deterministic central-exchange and bath-bath interactions
nb=4; sys.isotopes=repmat({'13C'},1,nb+1); sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,nb+1); inter.coupling.matrix=cell(nb+1);
spin_system=create(sys,inter); fields=[31 17 23 29 37]; tensors=cell(nb+1);
for j=1:nb+1
    spin_system.inter.zeeman.matrix{j}(3,3)=fields(j);
    for k=j+1:nb+1
        tensors{j,k}=[7 2 -3;-1 11 4;3 -2 13]/(j+k);
        spin_system.inter.coupling.matrix{j,k}=tensors{j,k}.';
    end
end
p.central_spins=1; p.bath_spins=2:5; p.method='gcce'; p.expansion='standard';
p.quantity='density'; p.order=4; v=[sqrt(0.3);1i*sqrt(0.7)]; p.central.rho=v*v';
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm'); p.meanfield='none';
p.times=linspace(0,0.07,8); p.sequence=[];
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^29);
maximum=0; trace_error=0; hermiticity=0; mineig=1;
for pure=[false true]
    states=repmat({eye(2)/2},1,nb);
    if pure
        for j=1:nb
            v=[cos(j/5);exp(1i*j/4)*sin(j/5)]; states{j}=v*v';
        end
    end
    p.bath_state=struct('mode','explicit_product','rho',{states});
    for hahn=[false true]
        p.sequence=[];
        if hahn
            p.sequence=struct('type','ideal','target',1,'placement','fraction',...
                              'at',0.5,'axis','y','angle',pi);
        end
        actual=cce(spin_system,p);
        for record=actual.raw{1}

            % Embed operators independently in descending physical spin order
            ids=[1 record.sites]; n=numel(ids); dim=2^n; db=dim/2;
            local={[0 1;1 0]/2,[0 -1i;1i 0]/2,diag([1 -1])/2};
            ops=cell(n,3); H=zeros(dim); rho=p.central.rho;
            for j=1:n
                for a=1:3
                    ops{j,a}=kron(eye(2^(j-1)),kron(local{a},eye(2^(n-j))));
                end
                H=H+fields(ids(j))*ops{j,3};
                if j>1, rho=kron(rho,states{ids(j)-1}); end
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
            pulse=expm(-1i*pi*ops{1,2}); reference=zeros(numel(p.times),4);
            for k=1:numel(p.times)
                t=p.times(k); u=expm(-1i*H*t);
                if hahn&&t>0, u=expm(-1i*H*t/2)*pulse*expm(-1i*H*t/2); end
                evolved=u*rho*u'; central=zeros(2);
                for a=1:2
                    for b=1:2
                        central(a,b)=trace(evolved((a-1)*db+(1:db),(b-1)*db+(1:db)));
                    end
                end
                reference(k,:)=central(:).';
            end
            discrepancy=max(abs(record.value-reference),[],'all'); maximum=max(maximum,discrepancy);
            result=test_true(result,'raw full-register exact',...
                discrepancy<=max(record.abs_error_bound,[],'all')+100*dim^3*eps,'Independent analytic or dense reference.');
        end
        trace_error=max(trace_error,max(actual.diagnostics.trace_error));
        hermiticity=max(hermiticity,max(actual.diagnostics.hermiticity_error));
        mineig=min(mineig,min(actual.diagnostics.minimum_eigenvalue));
    end
end
metrics.maximum_raw_error=maximum; metrics.trace_error=trace_error;
metrics.hermiticity_error=hermiticity; metrics.minimum_eigenvalue=mineig;

% Report the measured metrics
disp(metrics);

% Keep observable norm factors in output-column order on every evaluator path
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'processes',1};
inter=struct(); inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([2 1 5])/(2*pi);
inter.coupling.matrix{1,3}=diag([1 3 4])/(2*pi); small=create(sys,inter);
q=p; q.bath_spins=[2 3]; q.quantity='observable'; q.central.rho=diag([1 0]);
q.operators={eye(2),2*eye(2)}; q.bath_state=struct('mode','mixed_product');
q.jumps=struct('support',1,'operator',diag([1 -1]),'rate',0.7);
accepted=true; value_error=0; bound_error=0; cases=0;
finite=struct('type','finite','target',1,'placement','absolute','at',0.05,...
              'duration',0.05,'axis','x','amplitude',2,'detuning',0);
for method={'gcce','megcce'}
    q.method=method{1};
    for sequence={[],struct('name','hahn','count',1,'axis','x','target',1),finite}
        q.sequence=sequence{1};
        if strcmp(q.method,'megcce')&&~isempty(q.sequence), continue; end
        for mode={'mixed_product','interlaced'}
            q.bath_state=struct('mode',mode{1});
            if strcmp(mode{1},'interlaced')
                q.bath_state.samples=1; q.bath_state.seed=31;
                q.bath_state.probabilities='uniform'; q.bath_state.basis='uniform';
                q.interlaced.max_states=4;
            end
            for order=1:2
                q.order=order;
                for nt=2:3
                    q.times=linspace(0,0.5,nt);
                    for normalized=[false true]
                        q.normalized=normalized; q.operators={eye(2),2*eye(2)};
                        expected=column_evaluate(small,q); q.operators=q.operators.'; cases=cases+1;
                        try
                            actual=column_evaluate(small,q);
                            value_error=max(value_error,norm(actual.value-expected.value,'fro'));
                            bound_error=max(bound_error,max(abs(actual.abs_error_bound-...
                                            expected.abs_error_bound)./expected.abs_error_bound,[],'all'));
                        catch exception
                            accepted=false; fprintf('Column observable exception %s\n',exception.identifier);
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

% Mixed observable units must retain the same normalization decisions in column cells
observe=q; observe.operators={1e-20*eye(2),1e20*eye(2)}; observe.normalized=true;
row=cce_observe(q.central.rho,observe); observe.operators=observe.operators.';
column=NaN(size(row)); scaled_ok=true;
try
    column=cce_observe(q.central.rho,observe);
catch
    scaled_ok=false;
end
result=test_true(result,'Column scaled normalization accepted',scaled_ok,...
                 'Each reference must be compared only with its own operator scale.');
result=test_close(result,'Column scaled normalization values',column,row,1e-14,1e-14,...
                  'Independent unit changes commute with cell transposition.');

% An unresolved first observable must not be hidden by another valid column
observe.operators={diag([0 1]);eye(2)}; rejected=false;
try
    cce_observe(q.central.rho,observe);
catch exception
    rejected=strcmp(exception.identifier,'cce_observe:initial_zero');
end
result=test_true(result,'Column zero normalization rejected',rejected,...
                 'Every observable reference must be resolved before normalization.');
fprintf('R6 column-cell cases %d; accepted %d; value error %.3g; relative bound error %.3g\n',...
        cases,accepted,value_error,bound_error);

end

% Use the raw interlaced evaluator because its driver requires scalar coherence
function value=column_evaluate(spin_system,parameters)

% Exercise the public driver for supported non-interlaced observables
if ~strcmp(parameters.bath_state.mode,'interlaced')
    value=cce(spin_system,parameters); return;
end

% Evaluate the conditional raw average with matching scalar driver semantics
model=cce_prepare(spin_system,parameters); ensemble=cce_states(model,parameters);
schedule=cce_sequence(model,parameters);
parameters.interlaced.raw_only=true; parameters.interlaced.reference_state=1;
cluster=struct('id',uint64(1),'sites',parameters.bath_spins(1:parameters.order));
raw=cce_interlace(model,cluster,ensemble,schedule,parameters);
value=struct('value',raw.value,'abs_error_bound',raw.abs_error_bound);

end


