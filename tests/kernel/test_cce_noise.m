% Test: kernel/cce_noise | Additive bath noise
% Additive noise expansion against dense quantum and classical references. Syntax:
%
%                    result=test_cce_noise()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_noise()

% Prepare two bath spins and a nontracial product state
result=new_test_result('kernel/cce_noise','Additive bath noise',...
                       'Full order, disconnected cumulants, ordering and commuting limit.');
sys.isotopes={'1H','1H','1H'};
sys.magnet=0;
sys.output='hush';
inter.zeeman.scalar={0,0,0};
inter.coupling.scalar=cell(3);
s=create(sys,inter);
s.inter.zeeman.matrix{2}(3,3)=1.3;
s.inter.zeeman.matrix{3}(3,3)=0.7;
s=cce_test_write(s,[2 3],diag([0.4 0.5 0.9]));
p=struct('observed',1,'output','correlation','times',linspace(0,2,31),...
         'orientation',[0 0 0]);
p.cce=struct('bath_spins',[2 3],'expansion','standard','order',2,...
             'interaction_cutoff',Inf,'meanfield','none');
p.cce.graph=struct('mode','complete');
p.cce.bath_state=struct('mode','mixed_product');
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e8);
x=[0 1;1 0]/2;
y=[0 -1i;1i 0]/2;
z=diag([1 -1])/2;
p.noise=struct('ordering','ordered','connected',true,'manifold','bath_only',...
               'terms',struct('support',{2,3},'operator',{x,2*z}));
[r,diagnostics]=clexpn(s,@noise_correlation,p,'labframe');
B=kron(x,eye(2))+2*kron(eye(2),z);
H=1.3*kron(z,eye(2))+0.7*kron(eye(2),z)+...
  0.4*kron(x,x)+0.5*kron(y,y)+0.9*kron(z,z);
reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    U=expm(-1i*H*p.times(k));
    reference(k)=trace(U'*B*U*B)/4;
end
result=test_close(result,'C7 numerical check',r,reference,1e-12,1e-12,'Full-order additive exact reference.');
exact_error=max(abs(r-reference));

% Honour an explicitly supplied bath state independently of parameter defaults
model=cce_compile(s,p,'labframe');
state_parameters=diagnostics.parameters;
state_parameters.rho0=eye(2)/2;
ensemble=diagnostics.ensemble;
spin_config=rmfield(ensemble.spin_config(1),'assignment');
observed_state=diag([0 1]);
spin_config.bath_rho={diag([0.8 0.2]),[1 1;1 1]/2};
spin_config.expectations=[0 0 0.3;0.5 0 0];
contract=cce_capabilities(@noise_correlation);
schedule=struct('pulse_sequence',@noise_correlation,'contract',contract);
cluster=struct('id',1,'sites',[2 3]);
rho=kron(spin_config.bath_rho{1},spin_config.bath_rho{2});

% Direct covariance evaluation must bound retained records, not only operators
q=state_parameters;
q.rho0=observed_state;
q.times=linspace(0,1,128);
local=cce_subsystem(model,cluster.sites,spin_config,[]);
[ham,R,K,q]=cce_generators(local,q,[]);
q.cce.numerics.max_memory_bytes=4096;
refused=false;
try
    noise_correlation(local,q,ham,R,K);
catch exception
    refused=strcmp(exception.identifier,'noise_correlation:memory');
end
result=test_true(result,'Direct noise output memory preflight',refused,...
                 'The covariance and first-moment arrays must fit the declared byte budget.');
for ordering={'ordered','symmetrized'}
    q=state_parameters;
    q.rho0=observed_state;
    q.noise.ordering=ordering{1};
    q.noise.connected=true;
    raw=contract.evaluator(model,cluster,spin_config,schedule,q);
    supplied_reference=zeros(numel(p.times),1);
    for k=1:numel(p.times)
        U=expm(-1i*H*p.times(k));
        b_time=U'*B*U;
        greater=trace(b_time*B*rho);
        lesser=trace(B*b_time*rho);
        if strcmp(ordering{1},'ordered')
            moment=greater;
        else
            moment=(greater+lesser)/2;
        end
        supplied_reference(k)=moment-trace(b_time*rho)*trace(B*rho);
    end
    result=test_close(result,['Noise supplied spin_config ' ordering{1}],raw.value,supplied_reference,...
                      2e-12,2e-12,'Noise uses the spin_config bath state and its connected reference moments.');
end

% Check factorisation into independent connected single-spin correlations
s=cce_test_write(s,[2 3],zeros(3));
p.cce.order=1;
first=clexpn(s,@noise_correlation,p,'labframe');
p.cce.order=2;
second=clexpn(s,@noise_correlation,p,'labframe');
reference=cos(1.3*p.times(:))/4+1;
result=test_close(result,'C7 numerical check',first,reference,1e-12,1e-12,'Disconnected independent sum.');
result=test_close(result,'C7 numerical check',second,first,1e-12,1e-12,'Disconnected pair cumulant vanishes.');

% Exercise unequal-time quantum ordering in a polarised bath
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2]),eye(2)/2}});
p.noise.terms=struct('support',2,'operator',x);
ordered=clexpn(s,@noise_correlation,p,'labframe');
reference=(cos(1.3*p.times(:))+0.6i*sin(1.3*p.times(:)))/4;
result=test_close(result,'C7 numerical check',ordered,reference,1e-12,1e-12,'Ordered complex quantum correlation.');
p.noise.ordering='symmetrized';
symmetric=clexpn(s,@noise_correlation,p,'labframe');
result=test_close(result,'C7 numerical check',symmetric,real(reference),1e-12,1e-12,'Explicit anticommutator ordering.');

% Recover classical static variance for a commuting diagonal observable
p.noise.terms=struct('support',2,'operator',z);
classical=clexpn(s,@noise_correlation,p,'labframe');
result=test_close(result,'C7 numerical check',classical,repmat(0.16,numel(p.times),1),1e-12,1e-12,...
                  'Commuting classical connected variance.');

% Retain classical ensemble variance in a sampled connected correlation
q=p;
q.cce.bath_state=struct('mode','sampled','samples',64,'seed',31,...
                        'basis','uniform','probabilities','uniform');
q.rho0=eye(2)/2;
[~,diagnostics]=clexpn(s,@noise_correlation,q,'labframe');
ensemble=diagnostics.ensemble;
assignments=vertcat(ensemble.spin_config.assignment);
means=1.5-assignments(:,1);
expected=0.25-(ensemble.weights(:).'*means)^2;
for mode={'after_expansion','internal_then_expansion'}
    q.cce.averaging=mode{1};
    [sampled,diagnostics]=clexpn(s,@noise_correlation,q,'labframe');
    if strcmp(q.cce.averaging,'internal_then_expansion')
        result=test_true(result,'Internal connected-noise uncertainty unavailable',...
                         all(isnan(diagnostics.reduction.sampling_standard_error),'all'),...
                         'A connected sample covariance is not a linear sample mean with zero uncertainty.');
    end
    result=test_close(result,['Sampled connected noise ' mode{1}],sampled,...
                      repmat(expected,numel(q.times),1),1e-12,1e-12,...
                      'The mixed-state covariance includes the variance of sampled projected means.');
end

% Compare time-dependent ordered and symmetrised covariances with a dense mixture
s.inter.zeeman.matrix{2}=zeros(3);
s.inter.zeeman.matrix{2}(3,1)=1.3;
rho=diag([mean(assignments(:,1)==1) mean(assignments(:,1)==2)]);
for ordering={'ordered','symmetrized'}
    q.noise.ordering=ordering{1};
    sampled_dynamic=clexpn(s,@noise_correlation,q,'labframe');
    reference=zeros(size(sampled_dynamic));
    for k=1:numel(q.times)
        U=expm(-1i*1.3*x*q.times(k));
        bt=U'*z*U;
        moment=trace(bt*z*rho);
        if strcmp(ordering{1},'symmetrized')
            moment=trace((bt*z+z*bt)*rho)/2;
        end
        reference(k)=moment-trace(bt*rho)*trace(z*rho);
    end
    result=test_close(result,['Sampled dynamic covariance ' ordering{1}],sampled_dynamic,...
                      reference,1e-12,1e-12,'Time-dependent means must be averaged before connected subtraction.');
    fprintf('Dynamic %s covariance discrepancy %.3g.\n',ordering{1},max(abs(sampled_dynamic-reference)));
end
fprintf('R1236 sampled covariance %.16g; empirical variance %.16g.\n',...
        real(sampled(1)),expected);

fprintf('Noise full-order error %.3g; quantum ordering error %.3g\n',...
        exact_error,max(abs(ordered-reference)));

% Refuse routing a master-equation request through the noise evaluator
q=p;
q.rho0=ones(2)/2;
q.cce.basis=struct('formalism','zeeman-liouv');
clean=rmfield(q,'noise');
clean.output='signal';
model=cce_compile(s,clean,'labframe');
model.formalism='zeeman-liouv';
q.cce=model.parameters.cce;
contract=cce_capabilities(@noise_correlation);
schedule=struct('pulse_sequence',@free_evolution,'contract',cce_capabilities(@free_evolution));
cluster=struct('id',uint64(1),'sites',[2 3]);
foreign_config=struct('bath_ids',[2 3],'bath_rho',{p.cce.bath_state.rho},'weight',1);
rejected=false;
try
    contract.evaluator(model,cluster,foreign_config,schedule,q);
catch exception
    rejected=strcmp(exception.identifier,'clexpn:noise');
    fprintf('Foreign experiment validation: %s: %s\n',exception.identifier,exception.message);
end
result=test_true(result,'Noise evaluator rejects foreign method',rejected,...
                 'Direct evaluator calls must not silently use Hamiltonian noise dynamics for a master-equation method.');

% Reject non-Hermitian bath operators independently of their physical units
q=p;
q.cce.graph=struct('mode','complete');
q.noise.terms(1).operator=1e-20*[0 1;0 0];
rejected=false;
try
    clexpn(s,@noise_correlation,q,'labframe');
catch exception
    rejected=contains(lower(exception.message),'hermitian');
    fprintf('Scaled noise validation: %s\n',exception.message);
end
result=test_true(result,'Scaled non-Hermitian noise rejected',rejected,...
                 'A noise operator must remain Hermitian when its units are rescaled.');
fprintf('R6 scaled noise rejection %d\n',rejected);

% Retain disconnected observable cross terms in sparse graph families
s=create(sys,inter);
p.cce.order=2;
p.times=[0 0.3 0.8];
p.cce.graph=struct('mode','coupling','threshold',0);
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([1 0]),diag([1 0])}});
p.noise=struct('ordering','ordered','connected',false,'manifold','bath_only',...
               'terms',struct('support',{2,3},'operator',{z,z}));
uncoupled=clexpn(s,@noise_correlation,p,'labframe');
result=test_close(result,'disconnected unconnected noise',uncoupled,...
                  ones(numel(p.times),1),1e-12,0,...
                  'The two polarized singleton moments require their disconnected cross terms.');
fprintf('R7 disconnected noise zero-time %.17g\n',real(uncoupled(1)));

% Check evolving disconnected means against an independent dense reference
s.inter.zeeman.matrix{2}(3,1)=1.3;
s.inter.zeeman.matrix{3}(3,1)=0.7;
H=1.3*kron(x,eye(2))+0.7*kron(eye(2),x);
B=kron(z,eye(2))+kron(eye(2),z);
rho=kron(diag([1 0]),diag([1 0]));
reference=zeros(numel(p.times),1);
connected=reference;
for k=1:numel(p.times)
    U=expm(-1i*H*p.times(k));
    b_time=U'*B*U;
    reference(k)=trace(b_time*B*rho);
    connected(k)=reference(k)-trace(b_time*rho)*trace(B*rho);
end
for ordering={'ordered','symmetrized'}
    p.noise.ordering=ordering{1};
    p.noise.connected=false;
    actual=clexpn(s,@noise_correlation,p,'labframe');
    p.noise.connected=true;
    covariance=clexpn(s,@noise_correlation,p,'labframe');
    fprintf('R7 noise %s unconnected %.17g connected %.17g reference %.17g covariance %.17g\n',...
            ordering{1},actual(end),covariance(end),reference(end),connected(end));
    result=test_close(result,['sparse unconnected ' ordering{1}],actual,...
                      reference,1e-12,0,'Disconnected time-dependent means survive additive reduction.');
    result=test_close(result,['sparse connected ' ordering{1}],covariance,...
                      connected,1e-12,0,'Connected reduction still removes the disconnected moments.');
end

% Preserve cross-site covariance of the actual finite sampled mixture
p.cce.bath_state=struct('mode','sampled','samples',16,'seed',31,...
                        'basis','uniform','probabilities','uniform');
p.noise.connected=true;
p.rho0=eye(2)/2;
[~,diagnostics]=clexpn(s,@noise_correlation,p,'labframe');
ensemble=diagnostics.ensemble;
assignments=vertcat(ensemble.spin_config.assignment);
rho=zeros(4);
for k=1:size(assignments,1)
    index=2*(assignments(k,1)-1)+assignments(k,2);
    rho(index,index)=rho(index,index)+ensemble.weights(k);
end
for k=1:numel(p.times)
    U=expm(-1i*H*p.times(k));
    b_time=U'*B*U;
    reference(k)=trace(b_time*B*rho)-trace(b_time*rho)*trace(B*rho);
end
for averaging={'after_expansion','internal_then_expansion'}
    p.cce.averaging=averaging{1};
    actual=clexpn(s,@noise_correlation,p,'labframe');
    result=test_close(result,['sampled sparse covariance ' averaging{1}],...
                      actual,reference,1e-12,0,...
                      'The finite sampled mixture retains disconnected cross-site covariance.');
    fprintf('R7 sampled covariance discrepancy %.17g\n',max(abs(actual-reference)));
end

end


