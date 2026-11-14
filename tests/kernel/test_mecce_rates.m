% Analytic Lindblad rates, vectorization and conditional limitations.
%
% ilya.kuprov@weizmann.ac.il

function result=test_mecce_rates()

% Prepare a precessing central spin without a bath
result=new_test_result('kernel/mecce_rates','ME analytic rates',...
                       'Amplitude damping, dephasing, column convention and T1 limitation.');
sys.isotopes={'1H'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0}; s=create(sys,inter); s.inter.zeeman.matrix{1}(3,3)=2.3;
p=struct('central_spins',1,'bath_spins',[],'method','megcce','quantity','density',...
         'expansion','standard','order',1,'times',linspace(0,2,31),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,'meanfield','none');
p.graph=struct('mode','complete'); p.bath_state=struct('mode','mixed_product');
p.sequence=struct([]); p.central=struct('alpha',[1;0],'beta',[0;1],'rho',[1 1;1 1]/2);
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e8);
rate=0.7; p.jumps=struct('support',1,'operator',[0 1;0 0],'rate',rate);
r=cce(s,p); expected=exp(-(rate/2+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',squeeze(r.density(1,2,:)),expected(:),1e-12,1e-12,'Damping coherence rate and sign.');
result=test_close(result,'C7 numerical check',squeeze(r.density(2,2,:)),exp(-rate*p.times(:))/2,1e-12,1e-12,'Excited population decay.');

% Check the dephasing convention D[sigma_z] at the declared rate
p.jumps.operator=diag([1 -1]); r=cce(s,p);
expected=exp(-(2*rate+2.3i)*p.times)/2;
result=test_close(result,'C7 numerical check',squeeze(r.density(1,2,:)),expected(:),1e-12,1e-12,'Dephasing rate equals twice coefficient.');
p.method='mecce'; p.quantity='coherence'; conditional=cce(s,p);
result=test_close(result,'C7 numerical check',conditional.value,2*expected(:),1e-12,1e-12,'Conditional diagonal central jump.');

% Refuse nonclosed conditional jump blocks and negative rates
p.jumps.operator=[0 1;0 0]; m=cce_prepare(s,p); refused=false;
try
    cce_liouvillian(m,[],struct(),p);
catch err
    refused=strcmp(err.identifier,'cce:jumpBranches');
end
result=test_true(result,'C7 contract check',refused,'Central-changing jump requires ME-gCCE.');
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
                  'A vanishing coefficient contributes no dissipator or branch leakage.');

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
        max(abs(conditional.value-2*expected(:))),max(abs(first.*second-combined)));

end


