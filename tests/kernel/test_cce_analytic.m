% Analytic projected CCE limits, Hahn modulation, and complex FID. Syntax:
%
%                    result=test_cce_analytic()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_analytic()

% State independent targets rather than comparing two expansion routes
result=new_test_result('kernel/cce_analytic','CCE analytic limits',...
                       'Zero bath, isolated phase and one-spin Hahn modulation.');

% Build a spin-half pair with explicitly prescribed angular-frequency terms
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=zeros(3);
spin_system=create(sys,inter);
omega=1700;
longitudinal=700;
transverse=430;
phase=230;
spin_system.inter.zeeman.matrix{1}(3,3)=phase;
spin_system.inter.zeeman.matrix{2}(3,3)=omega;
spin_system.inter.coupling.matrix{1,2}(1,3)=transverse;
spin_system.inter.coupling.matrix{1,2}(3,3)=longitudinal;

% Declare the observed transition and bath policy
p.observed=1;
p.cce.bath_spins=2;
p.cce.approximation='projected';
p.cce.expansion='standard';
p.output='signal';
p.cce.order=1;
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','complete','score','coupling_norm');
p.cce.bath_state=struct('mode','mixed_product');
p.cce.meanfield='none';
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
p.times=linspace(0,0.012,81);

% Compare the complex projected FID against independent two-by-two expm
fid=clexpn(spin_system,@free_evolution,p,'labframe');
z=diag([1 -1])/2;
x=[0 1;1 0]/2;
ha=omega*z+(longitudinal*z+transverse*x)/2;
hb=omega*z-(longitudinal*z+transverse*x)/2;
reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    t=p.times(k);
    reference(k)=exp(-1i*phase*t)*trace(expm(-1i*ha*t)*expm(1i*hb*t))/2;
end
result=test_close(result,'complex FID',fid,reference,1e-11,1e-11,'Independent analytic or dense reference.');

% Detect one complex projected coherence with and without reference division
q=p;
samples=[1 4 9];
q.times=p.times(samples);
q.rho0=[0.5 -0.5i;0.5i 0.5];
q.coil=(2-3i)*[0 1;0 0];
initial=trace(q.coil'*q.rho0);
q.cce.numerics.atol=1e-30;
q.cce.numerics.rtol=1e-30;
frequency=norm(ha+(phase/2)*eye(2),1)+norm(hb-(phase/2)*eye(2),1);
for normalised=[false true]

    % Compare the signal and receiver-scaled roundoff against the one-spin formulae
    q.normalized=normalised;
    factor=initial;
    if normalised
        factor=1;
    end
    [signal,d]=clexpn(spin_system,@free_evolution,q,'labframe');
    model=cce_compile(spin_system,d.parameters,'labframe');
    check=cce_diagnostics(model,d.raw,d.reduction,d.parameters);
    expected=16*2^3*abs(factor)*((q.times(:)>0)+q.times(:)*frequency);
    result=test_close(result,'Projected explicit receiver FID',signal,...
                      factor*reference(samples),1e-12,1e-12,'Complex Hilbert-Schmidt detection.');
    result=test_close(result,'Projected receiver estimate in epsilon units',...
                      d.raw{1}(end).abs_error_bound/eps,expected,0,1e-12,...
                      'The mixed bath and trace detector have reciprocal Frobenius norms.');
    result=test_true(result,'Projected full bath remains valid above strict budget',...
                     all(check.valid)&&isequal(check.numerical_budget_met,[true;false;false])&&...
                     isequal(d.reduction.abs_error_bound,d.raw{1}(end).abs_error_bound),...
                     'Only zero elapsed time requires no propagation.');

end

% Evaluate the exact single-spin Hahn formula for two non-collinear fields
echo=clexpn(spin_system,@hahn_echo,p,'labframe');
a=[transverse/2 0 omega+longitudinal/2];
b=[-transverse/2 0 omega-longitudinal/2];
modulation=norm(cross(a,b))^2/(norm(a)^2*norm(b)^2);
reference=1-2*modulation*sin(norm(a)*p.times(:)/4).^2.*sin(norm(b)*p.times(:)/4).^2;
result=test_close(result,'one-spin Hahn formula',echo,reference,1e-11,1e-11,'Independent analytic or dense reference.');
result=test_true(result,'zero time is unity',echo(1)==1,'Independent analytic or dense reference.');

% Include absolute-time pulses only in samples that reach their timestamp
q=p;
q.times=[0 0.0005 0.001 0.0015];
q.refocusing=struct('spins',1,'operator',pi*x,'time',0.001);
[absolute,diagnostics]=clexpn(spin_system,@hahn_echo,q,'labframe');
schedule=diagnostics.raw{1}(end);
result=test_true(result,'Absolute event inclusion per sample',...
                 isequal(schedule.pulse_applied,[false false true true]),...
                 'Zero and shorter samples omit the pulse; equality and later samples include it.');
values=zeros(numel(q.times),1);
for k=1:numel(q.times)
    t=q.times(k);
    a=min(t,0.001);
    b=max(0,t-0.001);
    ua=expm(-1i*hb*b)*expm(-1i*ha*a);
    ub=expm(-1i*ha*b)*expm(-1i*hb*a);
    values(k)=exp(-1i*phase*(a-b))*trace(ua*ub')/2;
    if t>=q.refocusing.time
        values(k)=conj(values(k));
    end
end
result=test_close(result,'Absolute pulse fixed receiver reference',absolute,values,...
                  2e-12,2e-12,'Each requested duration uses only pulses at or before its endpoint.');
fprintf('Absolute pulse dense discrepancy %.3g.\n',max(abs(absolute-values)));

% Verify isolated-centre phase appears once in an empty bath
p.cce.bath_spins=[];
p.cce.order=0;
empty=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_close(result,'empty bath phase',empty,exp(-1i*phase*p.times(:)),1e-11,1e-11,'Independent analytic or dense reference.');

% Verify a disconnected unpolarised bath changes no observed dynamics
p.cce.bath_spins=2;
p.cce.order=1;
spin_system.inter.coupling.matrix{1,2}=zeros(3);
uncoupled=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_close(result,'uncoupled bath',uncoupled,empty,1e-11,1e-11,'Independent analytic or dense reference.');

% Evaluate the conditional bath trace at every requested time
values=zeros(numel(p.times),1);
for k=1:numel(p.times)
    values(k)=exp(-1i*phase*p.times(k))*trace(expm(-1i*ha*p.times(k))*expm(1i*hb*p.times(k)))/2;
end
metrics.fid_error=max(abs(fid-values));
metrics.hahn_error=max(abs(echo-reference));

% Check the same absolute boundary in both generalised propagation paths
probe_system=spin_system;
probe_system.inter.zeeman.matrix{1}(3,3)=phase;
probe_system.inter.coupling.matrix{1,2}=zeros(3);
probe_system.inter.coupling.matrix{1,2}(1,3)=transverse;
probe_system.inter.coupling.matrix{1,2}(3,3)=longitudinal;
probe=q;
probe.output='density';
probe.rho0=ones(2)/2;
probe.cce.approximation='generalized';
probe.cce.basis.formalism='zeeman-hilb';
full_hamiltonian=phase*kron(z,eye(2))+kron(eye(2),omega*z)+...
                 kron(z,longitudinal*z+transverse*x);
initial=kron(probe.rho0,eye(2)/2);
pulse=kron(expm(-1i*pi*x),eye(2));
expected_density=complex(zeros(numel(probe.times),4));
for index=1:numel(probe.times)
    time=probe.times(index);
    before=min(time,probe.refocusing.time);
    propagator=expm(-1i*before*full_hamiltonian);
    if time>0&&time>=probe.refocusing.time
        propagator=pulse*propagator;
    end
    propagator=expm(-1i*(time-before)*full_hamiltonian)*propagator;
    rho=propagator*initial*propagator';
    observed=complex(zeros(2));
    for row=1:2
        for col=1:2
            observed(row,col)=trace(rho((row-1)*2+(1:2),(col-1)*2+(1:2)));
        end
    end
    expected_density(index,:)=observed(:).';
end
for propagation={'eig','expm'}
    probe.cce.numerics.propagator=propagation{1};
    observed_density=clexpn(probe_system,@hahn_echo,probe,'labframe');
    result=test_close(result,['Absolute generalized ' propagation{1}],observed_density,...
                      expected_density,2e-12,2e-12,...
                      'Both propagators omit a later pulse and include the endpoint pulse.');
end

% Verify projected detection against the independently propagated full density
result=test_close(result,'Projected absolute full-density receiver',absolute,...
                  2*expected_density(:,3),2e-12,2e-12,...
                  'The fixed transition detector agrees with independent full-space propagation.');

% Report the measured metrics
disp(metrics);

end


