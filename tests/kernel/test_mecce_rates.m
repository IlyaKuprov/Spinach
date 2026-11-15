% Analytic Lindblad rates, vectorization and projected limitations.
%
% ilya.kuprov@weizmann.ac.il

function result=test_mecce_rates()

% Prepare a precessing observed spin without a bath
result=new_test_result('kernel/mecce_rates','ME analytic rates',...
                       'Amplitude damping, dephasing, column convention and T1 limitation.');
sys.isotopes={'1H'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0}; s=create(sys,inter); s.inter.zeeman.matrix{1}(3,3)=2.3;
p=struct('observed_spins',1,'bath_spins',[],'method','megcce','quantity','density',...
         'expansion','standard','order',1,'times',linspace(0,2,31),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,'meanfield','none');
p.graph=struct('mode','complete'); p.bath_state=struct('mode','mixed_product');
p.sequence=struct([]); p.rho0=[1 1;1 1]/2; p.transition=struct('alpha',[1;0],'beta',[0;1]);
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e8);
rate=0.7; p.jumps=struct('support',1,'operator',[0 1;0 0],'rate',rate);
r=cce(s,p); expected=exp(-(rate/2+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',squeeze(r.density(1,2,:)),expected(:),1e-12,1e-12,'Damping signal rate and sign.');
result=test_close(result,'C7 numerical check',squeeze(r.density(2,2,:)),exp(-rate*p.times(:))/2,1e-12,1e-12,'Excited population decay.');

% Resolve a weak dissipator below the stock generator clean-up threshold
weak=p; weak.times=[0 1]; weak.jumps.rate=1e-5; weak.numerics.propagator='stock';
weak.rho0=diag([1 0]); weak.jumps.operator=[0 0;1 0];
weak.numerics.atol=1e-8; weak.numerics.rtol=1e-8;
small=s; small.inter.zeeman.matrix{1}=zeros(3);
small.tols.small_matrix=1; small.tols.prop_chop=1e-4;
clean=cce(small,weak); population_error=abs(clean.density(1,1,end)-exp(-1e-5));
result=test_true(result,'Stock ME clean-up error covered',...
                 population_error<=clean.abs_error_bound(end,1),...
                 'The reported bound must cover weak dissipators below stock prop_chop.');
result=test_close(result,'Stock ME weak damping retained',clean.density(1,1,end),exp(-1e-5),...
                  1e-12,1e-12,'A CCE propagation budget cannot silently discard a physical jump.');

% Check the dephasing convention D[sigma_z] at the declared rate
p.jumps.operator=diag([1 -1]); r=cce(s,p);
expected=exp(-(2*rate+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',squeeze(r.density(1,2,:)),expected(:),1e-12,1e-12,'Dephasing rate equals twice coefficient.');
p.method='mecce'; p.quantity='signal'; projected=cce(s,p);
result=test_close(result,'C7 numerical check',projected.value,2*expected(:),1e-12,1e-12,'Projected diagonal observed jump.');

% Refuse nonclosed projected jump blocks and negative rates
p.jumps.operator=[0 1;0 0]; m=cce_prepare(s,p); refused=false;
try
    cce_liouvillian(m,[],struct(),p);
catch err
    refused=strcmp(err.identifier,'cce:manifold_mixing');
end
result=test_true(result,'C7 contract check',refused,'Observed-spin-changing jump requires ME-gCCE.');
% Reject manifold leakage independently of reciprocal operator-rate scaling
scaled=p; scaled.jumps.operator=1e-16*scaled.jumps.operator;
scaled.jumps.rate=1e32*scaled.jumps.rate; refused=false;
try
    cce_liouvillian(m,[],struct(),scaled);
catch err
    refused=strcmp(err.identifier,'cce:manifold_mixing');
end
result=test_true(result,'Scaled manifold-changing observed refusal',refused,...
                 'Reciprocal operator-rate scaling preserves the dissipator and its manifold leakage.');

% Continue the invalid-rate validation on the original operator
p.jumps.rate=-1; refused=false;
try
    cce_liouvillian(m,[],struct(),p);
catch
    refused=true;
end
result=test_true(result,'C7 contract check',refused,'Negative native dissipative rates are refused.');

% Require a zero-rate nonclosed jump to disappear rather than change the model
p.jumps.rate=0; [zero_rate,~]=cce_liouvillian(m,[],struct(),p);
p.jumps=struct('support',{},'operator',{},'rate',{});
[no_jump,~]=cce_liouvillian(m,[],struct(),p);
result=test_close(result,'Zero-rate limit',zero_rate,no_jump,0,0,...
                  'A vanishing coefficient contributes no dissipator or manifold leakage.');

% Refuse a dissipatively changing exterior field rather than freezing it silently
sys.isotopes={'1H','1H'}; inter.zeeman.scalar={0,0}; bath_system=create(sys,inter);
p.bath_spins=2; p.jumps=struct('support',2,'operator',[0 1;0 0],'rate',1);
model=cce_prepare(bath_system,p); refused=false;
try
    cce_liouvillian(model,[],struct('exterior_ids',2,'exterior_expectations',[0 0 0.5]),p);
catch err
    refused=strcmp(err.identifier,'cce:exteriorJump');
end
result=test_true(result,'Evolving exterior refusal',refused,...
                 'Active jumps cannot be silently omitted from a declared static exterior field.');

% Demonstrate why products of affine relaxing populations are not general T1
first=(1+exp(-0.4*p.times(:)))/2; second=(1+exp(-0.7*p.times(:)))/2;
combined=(1+exp(-1.1*p.times(:)))/2;
result=test_true(result,'C7 contract check',max(abs(first.*second-combined))>0.05,...
                 'Population products differ from the independent-rate relaxation solution.');
fprintf('ME dephasing maximum error %.3g; T1 product counterexample %.3g\n',...
        max(abs(projected.value-2*expected(:))),max(abs(first.*second-combined)));

% Refuse a nonnormal jump whose loss operator leaves a selected manifold
sys.isotopes={'14N'}; inter.zeeman.scalar={0}; s=create(sys,inter);
a=[1;0;0]; b=[0;1;0]; psi=(a+b)/sqrt(2);
p.method='mecce'; p.quantity='signal'; p.times=[0 1]; p.bath_spins=[];
p.rho0=psi*psi'; p.transition=struct('alpha',a,'beta',b);
op=[1 0 1;0 0 0;0 0 0];
p.jumps=struct('support',1,'operator',op,'rate',1);
rejected=false;
try
    cce(s,p);
catch exception
    rejected=strcmp(exception.identifier,'cce:manifold_mixing');
end
result=test_true(result,'Projected loss-operator leakage refused',rejected,...
                 'Jump invariance alone does not close a nonnormal Lindblad signal block.');

% Retain the excluded level in an independent dense master-equation reference
p.method='megcce'; full=cce(s,p); loss=op'*op;
generator=kron(conj(op),op)-(kron(eye(3),loss)+kron(loss.',eye(3)))/2;
rho=reshape(expm(generator)*p.rho0(:),3,3);
result=test_close(result,'Nonnormal jump full-space reference',full.value(end),rho(1,2),...
                  2e-12,2e-12,'ME-gCCE retains the third-level loss-operator coupling.');
fprintf('R5567 normalized full %.16g; projected truncation %.16g.\n',...
        2*real(rho(1,2)),exp(-0.5));

% Continue accepting reducing projected subspaces with an excluded level
p.method='mecce'; p.jumps.operator=diag([1 0 2]); closed=cce(s,p);
result=test_close(result,'Reducing three-level projected jump',closed.value(end),exp(-0.5),...
                  2e-12,2e-12,'An excluded level alone is not grounds for refusing a closed block.');

end


