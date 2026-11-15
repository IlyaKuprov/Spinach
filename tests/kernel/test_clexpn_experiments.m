% Validate closed-system experiments against independent matrix propagation.
% Syntax:
%
%                      result=test_clexpn_experiments()
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_experiments()

% Noncommuting bath terms produce a complex, nontrivial echo modulation
result=new_test_result('kernel/clexpn_experiments','Native cluster experiments',...
                       'Free evolution, Hahn, CPMG, ESEEM and spectral bath pumping.');
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush'; sys.disable={'hygiene'};
inter.zeeman.matrix={diag([0 0 3]),diag([0 0 11])};
inter.coupling.matrix=cell(2); inter.coupling.matrix{1,2}=[0 0 0;0 0 0;5 0 7];
spin_system=create(sys,inter); p.observed=1; p.times=[0 0.007 0.019];
p.transition.spins=1; p.transition.alpha=[1;0]; p.transition.beta=[0;1];
p.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8);
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2])}});
[p,~]=cce_options(spin_system,p,@free_evolution,'labframe');
model=cce_compile(spin_system,p,'labframe');
[cluster,~]=cce_subsystem(model,2,struct(),[]);
[H,R,K,p]=cce_generators(cluster,p,[]);
[rho0,coil,reference]=cce_observation(cluster,p);
pulse=kron([0 -1i;-1i 0],eye(2));
sequences={@free_evolution,@hahn_echo,@cpmg_echo,@cce_eseem}; p.nloops=3;
fractions={[],0.5,([1 2 3]-0.5)/3,0.5};
for j=1:numel(sequences)
    actual=sequences{j}(cluster,p,H,R,K); expected=zeros(size(actual));
    for n=1:numel(p.times)
        u=eye(4); previous=0;
        if p.times(n)>0
            for fraction=fractions{j}
                u=pulse*expm(-1i*H*(fraction-previous)*p.times(n))*u; previous=fraction;
            end
        end
        u=expm(-1i*H*(1-previous)*p.times(n))*u;
        expected(n)=trace(coil'*u*rho0*u')/reference;
    end
    result=test_close(result,['Independent ' func2str(sequences{j})],actual,expected,3e-13,3e-13,...
                      'Dense expm independently checks elapsed time, pulses and complex receiver.');
end

% Empty spectral selection is a documented no-pump control
p.pump=struct('center',1e12,'halfwidth',0,'units','Hz','fraction',0.25,'allow_empty',true);
actual=bath_pumped_hahn(cluster,p,H,R,K); expected=hahn_echo(cluster,p,H,R,K);
result=test_close(result,'Empty allowed pump window',actual,expected,3e-13,3e-13,...
                  'An empty explicit window leaves the Hahn experiment unchanged.');

% Free-evolution density has the same column-major ordering as direct partial trace
p.output='density'; actual=free_evolution(cluster,p,H,R,K); expected=zeros(3,4);
for n=1:numel(p.times)
    u=expm(-1i*H*p.times(n)); rho=u*rho0*u'; observed=zeros(2);
    for a=1:2
        for b=1:2, observed(a,b)=trace(rho((a-1)*2+(1:2),(b-1)*2+(1:2))); end
    end
    expected(n,:)=observed(:).';
end
result=test_close(result,'Column-major reduced density',actual,expected,3e-13,3e-13,...
                  'Direct free evolution returns physical observed density elements.');

% A physical receiver phase remains complex before normalization
phase=0.37; phased=p; phased.output='observable'; phased.normalized=false;
phased.coil=exp(1i*phase)*model.transition.alpha*model.transition.beta';
unphased=phased; unphased.coil=model.transition.alpha*model.transition.beta';
actual=free_evolution(cluster,phased,H,R,K);
expected=exp(-1i*phase)*free_evolution(cluster,unphased,H,R,K);
result=test_close(result,'Complex receiver phase',actual,expected,3e-13,3e-13,...
                  'Detection uses the conjugate-transposed receiver operator.');

% Inclusive spectral-window boundaries select the same physical bath pulse
p.output='signal'; p.pump.center=model.zeeman{2}(3,3); p.pump.units='rad/s';
p.pump.halfwidth=0; p.pump.allow_empty=false;
actual=bath_pumped_hahn(cluster,p,H,R,K);
reference_pump=p; reference_pump.pump.halfwidth=eps(abs(p.pump.center));
expected=bath_pumped_hahn(cluster,reference_pump,H,R,K);
result=test_close(result,'Inclusive pump boundary',actual,expected,0,0,...
                  'A site exactly on the spectral boundary is retained.');

end


