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

end


