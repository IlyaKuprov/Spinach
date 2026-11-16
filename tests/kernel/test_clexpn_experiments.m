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

% Finite drives remain experiment-owned and retain the coupled drift
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'processes',1};
inter.zeeman.matrix={zeros(3),zeros(3)}; inter.coupling.matrix=cell(2);
system=create(sys,inter); p=struct(); p.observed=1; p.times=[0 0.2];
p.transition=struct('spins',1,'ms',[0.5 -0.5]); p.output='density';
p.cce=struct('approximation','generalized','order',1,'graph',struct('mode','complete'));
p.sequence=struct('type','finite','target',1,'axis','z',...
    'placement','absolute','at',0,'duration',0.2,'amplitude',7);
actual=clexpn(system,@free_evolution,p,'labframe');
rho=[1 1;1 1]/2; U=expm(-1i*7*diag([0.5 -0.5])*0.2); expected=U*rho*U';
result=test_close(result,'Finite-drive acceptance probe',actual(2,:),expected(:).',...
                  2e-14,0,'The supplied finite event rotates the observed density.');
p.preparation='transverse'; named=clexpn(system,@free_evolution,p,'labframe');
result=test_close(result,'Transverse preparation',named,actual,2e-14,0,...
                  'The named positive preparation matches the transition superposition.');
p.cce.basis.formalism='zeeman-liouv';
liouville=clexpn(system,@free_evolution,p,'labframe');
result=test_close(result,'Finite-drive Liouville representation',liouville,actual,...
                  2e-14,0,'The Hamiltonian commutator gives the same observed density.');
p.cce=rmfield(p.cce,'basis'); p.times=[0 0.05 0.2 0.4];
p.sequence.duration=0.1; p.sequence.at=0.05; p.sequence.detuning=2;
p.sequence(2)=p.sequence(1); p.sequence(2).at=0.1;
p.sequence(2).axis='x'; p.sequence(2).amplitude=3;
actual=clexpn(system,@free_evolution,p,'labframe'); expected=zeros(size(actual));
Sz=diag([0.5 -0.5]); Sx=[0 1;1 0]/2;
for n=1:numel(p.times)
    t=p.times(n); state=rho; cuts=unique([0 t min([0.05 0.1 0.15 0.2],t)]);
    for k=1:numel(cuts)-1
        middle=(cuts(k)+cuts(k+1))/2; generator=zeros(2);
        if middle>0.05&&middle<0.15, generator=generator+9*Sz; end
        if middle>0.1&&middle<0.2, generator=generator+3*Sx+2*Sz; end
        U=expm(-1i*generator*(cuts(k+1)-cuts(k))); state=U*state*U';
    end
    expected(n,:)=state(:).';
end
result=test_close(result,'Overlapping finite intervals',actual,expected,3e-14,0,...
                  'Start/stop boundaries, detuning and partial intervals match dense expm.');
p.sequence=p.sequence(1); p.sequence.axis=[0 0 1];
p.sequence.placement='fraction'; p.sequence.at=0.5;
actual=clexpn(system,@free_evolution,p,'labframe');
for n=1:numel(p.times)
    U=expm(-1i*9*Sz*min(0.1,p.times(n)/2)); state=U*rho*U';
    expected(n,:)=state(:).';
end
result=test_close(result,'Fractional finite-drive placement',actual,expected,3e-14,0,...
                  'The duration remains in seconds while start scales with elapsed time.');
% Noncommuting drift remains active throughout a finite control interval
inter.coupling.matrix{1,2}=diag([1 2 3]); system=create(sys,inter);
p.sequence.at=0; p.sequence.duration=0.1; p.sequence.axis='x';
p.sequence.placement='absolute'; p.preparation='transverse';
actual=clexpn(system,@free_evolution,p,'labframe');
Sy=[0 -1i;1i 0]/2;
drift=2*pi*(kron(Sx,Sx)+2*kron(Sy,Sy)+3*kron(Sz,Sz));
drive=kron(7*Sx+2*Sz,eye(2));
for n=1:numel(p.times)
    duration=min(0.1,p.times(n));
    U=expm(-1i*drift*(p.times(n)-duration))*expm(-1i*(drift+drive)*duration);
    state=U*kron(rho,eye(2)/2)*U'; observed=zeros(2);
    for a=1:2
        for b=1:2, observed(a,b)=trace(state((a-1)*2+(1:2),(b-1)*2+(1:2))); end
    end
    expected(n,:)=observed(:).';
end
result=test_close(result,'Finite-drive coupled drift',actual,expected,6e-14,0,...
                  'Dense full-bath propagation includes noncommuting dipolar drift.');

% Relaxation remains simultaneous with a finite observed drive
inter.coupling.matrix=cell(2); system=create(sys,inter);
p.cce.basis.formalism='zeeman-liouv';
p.bath_jumps=struct('support',1,'operator','Lz','rate',2);
actual=clexpn(system,@free_evolution,p,'labframe');
J=sqrt(2)*Sz; loss=J'*J;
D=kron(conj(J),J)-(kron(eye(2),loss)+kron(loss.',eye(2)))/2;
A=-1i*(kron(eye(2),7*Sx+2*Sz)-kron((7*Sx+2*Sz).',eye(2)))+D;
for n=1:numel(p.times)
    duration=min(0.1,p.times(n));
    expected(n,:)=(expm(D*(p.times(n)-duration))*expm(A*duration)*rho(:)).';
end
result=test_close(result,'Finite-drive GKSL evolution',actual,expected,6e-14,0,...
                  'The dissipator remains active during and after the drive.');
p=rmfield(p,'bath_jumps'); p.cce=rmfield(p.cce,'basis');

p.sequence.type='ideal';
failed=false;
try
    clexpn(system,@free_evolution,p,'labframe');
catch err
    failed=strcmp(err.identifier,'free_evolution:sequence');
end
result=test_true(result,'Free evolution rejects ideal events',failed,...
                 'Named echo experiments own ideal refocusing pulses.');
p.sequence.type='finite'; p.cce.approximation='projected'; p.output='signal';
failed=false;
try
    clexpn(system,@free_evolution,p,'labframe');
catch err
    failed=strcmp(err.identifier,'free_evolution:projection');
end
result=test_true(result,'Projected finite drives rejected',failed,...
                 'A finite observed drive requires generalized propagation.');


% Product-ms selectors prepare all high-spin and multiple-centre examples
isotopes={{'E4'},{'E3'},{'E','E'},{'E','51V'}};
selectors={[-0.5 0.5],[-1 1],[-0.5 0.5;-0.5 0.5],[-0.5 0.5;3.5 3.5]};
positions={[2 3],[1 3],[1 4],[1 9]}; dimensions=[4 3 4 16];
for k=1:numel(isotopes)
    sys.isotopes=isotopes{k}; inter=struct();
    inter.zeeman.matrix=repmat({zeros(3)},1,numel(sys.isotopes));
    inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes)); system=create(sys,inter);
    p=struct(); p.observed=1:numel(sys.isotopes); p.times=0; p.output='density';
    p.transition=struct('spins',p.observed,'ms',selectors{k});
    p.cce=struct('approximation','generalized','order',0);
    actual=clexpn(system,@free_evolution,p,'labframe');
    vector=zeros(dimensions(k),1); vector(positions{k})=1/sqrt(2);
    expected=vector*vector';
    result=test_close(result,sprintf('Example preparation selector %d',k),...
                      actual,expected(:).',2e-14,0,...
                      'Selected product manifolds replace hand-built preparation vectors.');
end


end


