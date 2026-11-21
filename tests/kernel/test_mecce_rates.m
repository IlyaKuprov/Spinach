% Test: kernel/mecce_rates | ME analytic rates
% Analytic Lindblad rates, vectorisation, and projected limitations.
% Density answers are time-by-column-vectorised observed matrices.
% Projected no-bath generators are scalars; an empty dissipator is zero.
%
% Syntax:
%
%                    result=test_mecce_rates()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_mecce_rates()

% Prepare a precessing observed spin without a bath
result=new_test_result('kernel/mecce_rates','ME analytic rates',...
                       'Amplitude damping, dephasing, column convention and T1 limitation.');
sys.isotopes={'1H'};
sys.magnet=0;
sys.output='hush';
inter.zeeman.scalar={0};
s=create(sys,inter);
s.inter.zeeman.matrix{1}(3,3)=2.3;
p.observed=1;
p.output='density';
p.times=linspace(0,2,31);
p.orientation=[0 0 0];
p.cce=struct('bath_spins',[],'approximation','generalized',...
             'expansion','standard','order',0,'interaction_cutoff',Inf,'meanfield','none');
p.cce.graph=struct('mode','complete');
p.cce.bath_state=struct('mode','mixed_product');
p.cce.basis=struct('formalism','zeeman-liouv');
p.rho0=[1 1;1 1]/2;
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e8);
rate=0.7;
p.bath_jumps=struct('support',1,'operator',[0 1;0 0],'rate',rate);
r=clexpn(s,@free_evolution,p,'labframe');
expected=exp(-(rate/2+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',r(:,3),expected(:),1e-12,1e-12,'Damping signal rate and sign.');
result=test_close(result,'C7 numerical check',r(:,4),exp(-rate*p.times(:))/2,1e-12,1e-12,'Excited population decay.');

% Resolve a weak dissipator below the stock generator clean-up threshold
weak=p;
weak.times=[0 1];
weak.bath_jumps.rate=1e-5;
weak.cce.numerics.propagator='stock';
weak.rho0=diag([1 0]);
weak.bath_jumps.operator=[0 0;1 0];
weak.cce.numerics.atol=1e-8;
weak.cce.numerics.rtol=1e-8;
small=s;
small.inter.zeeman.matrix{1}=zeros(3);
small.bas.formalism='zeeman-liouv';
small.tols.small_matrix=1;
small.tols.prop_chop=1e-4;
weak_model=cce_compile(small,weak,'labframe');
[weak_system,~]=cce_subsystem(weak_model,[],struct(),[]);
[H,R,K,weak]=cce_generators(weak_system,weak_model.parameters,[]);
[clean,clean_details]=free_evolution(weak_system,weak,H,R,K);
population_error=abs(clean(end,1)-exp(-1e-5));
result=test_true(result,'Stock ME clean-up error covered',...
                 population_error<=clean_details.abs_error_bound(end,1),...
                 'The reported bound must cover weak dissipators below stock prop_chop.');
result=test_close(result,'Stock ME weak damping retained',clean(end,1),exp(-1e-5),...
                  1e-12,1e-12,'A CCE propagation budget cannot silently discard a physical jump.');

% Check the dephasing convention D[2Sz] at the declared rate
p.bath_jumps.operator=diag([1 -1]);
r=clexpn(s,@free_evolution,p,'labframe');
expected=exp(-(2*rate+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',r(:,3),expected(:),1e-12,1e-12,...
                  'Dephasing rate equals twice coefficient.');
p.cce.approximation='projected';
p.output='signal';
projected=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'C7 numerical check',projected,2*expected(:),1e-12,1e-12,'Projected diagonal observed jump.');

% Refuse nonclosed projected jump blocks and negative rates
p.bath_jumps.operator=[0 1;0 0];
s.bas.formalism='zeeman-liouv';
m=cce_compile(s,p,'labframe');
p=m.parameters;
[subsystem,~]=cce_subsystem(m,[],struct(),[]);
refused=false;
try
    cce_generators(subsystem,p,[]);
catch err
    refused=strcmp(err.identifier,'cce:capability');
end
result=test_true(result,'C7 contract check',refused,'Observed-spin-changing jump requires ME-gCCE.');

% Reject manifold leakage independently of reciprocal operator-rate scaling
scaled=p;
scaled.bath_jumps.operator=1e-16*scaled.bath_jumps.operator;
scaled.bath_jumps.rate=1e32*scaled.bath_jumps.rate;
refused=false;
try
    cce_generators(subsystem,scaled,[]);
catch err
    refused=strcmp(err.identifier,'cce:capability');
end
result=test_true(result,'Scaled manifold-changing observed refusal',refused,...
                 'Reciprocal operator-rate scaling preserves the dissipator and its manifold leakage.');

% Continue the invalid-rate validation on the original operator
p.bath_jumps.rate=-1;
refused=false;
try
    cce_generators(subsystem,p,[]);
catch
    refused=true;
end
result=test_true(result,'C7 contract check',refused,'Negative native dissipative rates are refused.');

% Require a zero-rate nonclosed jump to disappear rather than change the model
p.bath_jumps.rate=0;
zero_rate=cce_generators(subsystem,p,[]);
p.bath_jumps=struct('support',{},'operator',{},'rate',{});
no_jump=cce_generators(subsystem,p,[]);
zero_loss=0;
no_loss=0;
if ~isempty(zero_rate.dissipator)
    zero_loss=zero_rate.dissipator;
end
if ~isempty(no_jump.dissipator)
    no_loss=no_jump.dissipator;
end
zero_rate=-1i*(zero_rate.left-zero_rate.right)+zero_loss;
no_jump=-1i*(no_jump.left-no_jump.right)+no_loss;
result=test_close(result,'Zero-rate limit',zero_rate,no_jump,0,0,...
                  'A vanishing coefficient contributes no dissipator or manifold leakage.');

% Refuse a dissipatively changing exterior field rather than freezing it silently
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
bath_system=create(sys,inter);
p.cce.bath_spins=2;
p.bath_jumps=struct('support',2,'operator',[0 1;0 0],'rate',1);
p.cce.order=1;
model=cce_compile(bath_system,p,'labframe');
model.parameters.cce.meanfield='static';
refused=false;
try
    configuration=struct('bath_ids',2,'bath_rho',{{diag([1 0])}},...
                         'exterior_ids',2,'exterior_expectations',[0 0 0.5]);
    [exterior_system,~]=cce_subsystem(model,[],configuration,[]);
    cce_generators(exterior_system,model.parameters,[]);
catch err
    refused=strcmp(err.identifier,'cce:exteriorJump');
end
result=test_true(result,'Evolving exterior refusal',refused,...
                 'Active jumps cannot be silently omitted from a declared static exterior field.');

% Demonstrate why products of affine relaxing populations are not general T1
first=(1+exp(-0.4*p.times(:)))/2;
second=(1+exp(-0.7*p.times(:)))/2;
combined=(1+exp(-1.1*p.times(:)))/2;
result=test_true(result,'C7 contract check',max(abs(first.*second-combined))>0.05,...
                 'Population products differ from the independent-rate relaxation solution.');
fprintf('ME dephasing maximum error %.3g; T1 product counterexample %.3g\n',...
        max(abs(projected-2*expected(:))),max(abs(first.*second-combined)));

% Refuse a nonnormal jump whose loss operator leaves a selected manifold
sys.isotopes={'14N'};
inter.zeeman.scalar={0};
s=create(sys,inter);
a=[1;0;0];
b=[0;1;0];
psi=(a+b)/sqrt(2);
p.cce.approximation='projected';
p.output='signal';
p.times=[0 1];
p.cce.bath_spins=[];
p.cce.order=0;
p.rho0=psi*psi';
p.transition=struct('spins',1,'alpha',a,'beta',b);
op=[1 0 1;0 0 0;0 0 0];
p.bath_jumps=struct('support',1,'operator',op,'rate',1);
rejected=false;
try
    clexpn(s,@free_evolution,p,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce:capability');
end
result=test_true(result,'Projected loss-operator leakage refused',rejected,...
                 'Jump invariance alone does not close a nonnormal Lindblad signal block.');

% Retain the excluded level in an independent dense master-equation reference
p.cce.approximation='generalized';
p.normalized=false;
full=clexpn(s,@free_evolution,p,'labframe');
loss=op'*op;
generator=kron(conj(op),op)-(kron(eye(3),loss)+kron(loss.',eye(3)))/2;
rho=reshape(expm(generator)*p.rho0(:),3,3);
result=test_close(result,'Nonnormal jump full-space reference',full(end),rho(1,2),...
                  2e-12,2e-12,'ME-gCCE retains the third-level loss-operator coupling.');
fprintf('R5567 normalised full %.16g; projected truncation %.16g.\n',...
        2*real(rho(1,2)),exp(-0.5));

% Continue accepting reducing projected subspaces with an excluded level
p.cce.approximation='projected';
p.normalized=true;
p.bath_jumps.operator=diag([1 0 2]);
closed=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Reducing three-level projected jump',closed(end),exp(-0.5),...
                  2e-12,2e-12,'An excluded level alone is not grounds for refusing a closed block.');

end


