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
p.drive=struct('spins',1,'axis','z','time',0,...
               'duration',0.2,'amplitude',7,'detuning',0);
actual=clexpn(system,@finite_drive,p,'labframe');
rho=[1 1;1 1]/2; U=expm(-1i*7*diag([0.5 -0.5])*0.2); expected=U*rho*U';
result=test_close(result,'Finite-drive acceptance probe',actual(2,:),expected(:).',...
                  2e-14,0,'The native finite drive rotates the observed density.');
p.preparation='transverse'; named=clexpn(system,@finite_drive,p,'labframe');
result=test_close(result,'Transverse preparation',named,actual,2e-14,0,...
                  'The named positive preparation matches the transition superposition.');
p.cce.basis.formalism='zeeman-liouv';
liouville=clexpn(system,@finite_drive,p,'labframe');
result=test_close(result,'Finite-drive Liouville representation',liouville,actual,...
                  2e-14,0,'The Hamiltonian commutator gives the same observed density.');
p.cce=rmfield(p.cce,'basis'); p.times=[0 0.05 0.2 0.4];
p.drive.duration=0.1; p.drive.time=0.05; p.drive.detuning=2;
p.drive(2)=p.drive(1); p.drive(2).time=0.1;
p.drive(2).axis='x'; p.drive(2).amplitude=3;
actual=clexpn(system,@finite_drive,p,'labframe'); expected=zeros(size(actual));
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
p.drive=p.drive(1); p.drive.axis=[0 0 1];
p.drive=rmfield(p.drive,'time'); p.drive.fraction=0.5;
actual=clexpn(system,@finite_drive,p,'labframe');
for n=1:numel(p.times)
    U=expm(-1i*9*Sz*min(0.1,p.times(n)/2)); state=U*rho*U';
    expected(n,:)=state(:).';
end
result=test_close(result,'Fractional finite-drive placement',actual,expected,3e-14,0,...
                  'The duration remains in seconds while start scales with elapsed time.');

% Noncommuting drift remains active throughout a finite control interval
inter.coupling.matrix{1,2}=diag([1 2 3]); system=create(sys,inter);
p.drive=rmfield(p.drive,'fraction');
p.drive.time=0; p.drive.duration=0.1; p.drive.axis='x';
p.preparation='transverse';
actual=clexpn(system,@finite_drive,p,'labframe');
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

% Named finite drives compile the same control without a public event array
q=rmfield(p,'drive');
q.drive=struct('spins',1,'axis','x','amplitude',7,'detuning',2,'duration',0.1);
selected=p.times==0|p.times>=q.drive.duration; q.times=p.times(selected);
named=clexpn(system,@finite_drive,q,'labframe');
result=test_close(result,'Named finite-drive contract',named,actual(selected,:),0,0,...
                  'The declarative drive preserves the controlled density bitwise.');
invalid={struct('spins',2),struct('axis',[1 1 0]),struct('amplitude',Inf),...
         struct('detuning',1i),struct('duration',0),struct('extra',1)};
for n=1:numel(invalid)
    bad=q; field=fieldnames(invalid{n}); bad.drive.(field{1})=invalid{n}.(field{1});
    rejected=false;
    try
        clexpn(system,@finite_drive,bad,'labframe');
    catch
        rejected=true;
    end
    result=test_true(result,['Named finite-drive validation ' field{1}],rejected,...
                      'Every declarative drive field is validated.');
end
bad=q; bad.times=[0 0.02]; partial=clexpn(system,@finite_drive,bad,'labframe');
U=expm(-1i*(drift+drive)*0.02); state=U*kron(rho,eye(2)/2)*U'; observed=zeros(2);
for a=1:2
    for b=1:2, observed(a,b)=trace(state((a-1)*2+(1:2),(b-1)*2+(1:2))); end
end
result=test_close(result,'Named finite-drive observation time',partial(2,:),observed(:).',...
                  6e-14,0,'An observation during a drive retains simultaneous coupled evolution.');
bad=q; bad.sequence=struct('type','finite','target',1,'axis','x',...
    'placement','absolute','at',0,'duration',0.1,'amplitude',7,'detuning',2); rejected=false;
try
    clexpn(system,@finite_drive,bad,'labframe');
catch
    rejected=true;
end
result=test_true(result,'Named finite-drive event ownership',rejected,...
                  'A competing event array is rejected.');

% Mixed finite and ideal controls retain their explicit representation limit
bad=q; bad.cce.basis.formalism='zeeman-liouv';
bad.refocusing={struct('spins',1,'axis','x','angle',pi,'fraction',0.5)};
rejected=false;
try
    clexpn(system,@finite_drive,bad,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_acquire:controls');
end
result=test_true(result,'Mixed Liouville controls rejected',rejected,...
                 'Ideal refocusing mixed with finite drives requires generalized Hilbert dynamics.');

% Relaxation remains simultaneous with a finite observed drive
inter.coupling.matrix=cell(2); system=create(sys,inter);
p.cce.basis.formalism='zeeman-liouv';
p.bath_jumps=struct('support',1,'operator','Lz','rate',2);
actual=clexpn(system,@finite_drive,p,'labframe');
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

% Raw finite, ideal, and empty event arrays are all removed inputs
for events={struct('type','finite'),struct('type','ideal'),struct([]),[]}
    rejected=p; rejected=rmfield(rejected,'drive'); rejected.sequence=events{1};
    failed=false;
    try
        clexpn(system,@free_evolution,rejected,'labframe');
    catch err
        failed=strcmp(err.identifier,'cce_options:sequence');
    end
    result=test_true(result,'Free evolution rejects raw events',failed,...
                     'Use finite_drive, refocusing_train, or a named echo experiment.');
end
p.cce.approximation='projected'; p.output='signal';
failed=false;
try
    clexpn(system,@finite_drive,p,'labframe');
catch err
    failed=strcmp(err.identifier,'clexpn:projection');
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


% Coherent diagnostics retain the original arithmetic and exact-reference masks
sys.magnet=0; sys.isotopes={'E','1H'}; sys.output='hush';
inter=struct(); inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2); inter.coupling.matrix{1,2}=diag([1 2 3]);
system=create(sys,inter); p=struct();
p.observed=1; p.transition.spins=1; p.transition.ms=[0.5 -0.5]; p.times=[0 0.1];
p.output='density'; p.cce.approximation='generalized';
p.cce.numerics=struct('atol',1e-20,'rtol',1e-20,'max_memory_bytes',2^28);
[~,d]=clexpn(system,@free_evolution,p,'labframe');
model=cce_compile(system,d.parameters,'labframe');
[stock_transition,stock_meta]=cce_transition(system,d.parameters,'labframe');
result=test_close(result,'Compiled observed transition matches stock',...
                  model.transition.alpha,stock_transition.alpha,0,0,...
                  'Native observed-only preparation uses the same ordered physical state.');
result=test_close(result,'Compiled observed energies match stock',...
                  model.transition_meta.energies,stock_meta.energies,1e-13,1e-13,...
                  'Independent stock Hamiltonian assembly validates the observed generator.');
configuration=d.ensemble.spin_config; configuration=rmfield(configuration,'assignment');
[retained,~]=cce_subsystem(model,2,configuration,[]);
[H,R,K,pseq]=cce_generators(retained,d.parameters,[]);
[~,aux]=free_evolution(retained,pseq,H,R,K);
[v,e]=eig(full((H+H')/2),'vector');
residual=norm(H*v-v.*e.','fro'); orthogonality=norm(v'*v-eye(size(H)),'fro');
operations=12*size(H,1)^3; gamma=operations*eps/(1-operations*eps);
bound=2*p.times(:)*(residual+pseq.cce_aux.roundoff_bound)+...
      4*orthogonality+gamma*(1+p.times(:)*norm(H,'fro'));
result=test_close(result,'Native coherent error bounds',aux.abs_error_bound,...
                  repmat(bound,1,4),0,0,'The original propagation-error estimate is retained bitwise.');
result=test_close(result,'Measured spectral residual',aux.eigen_residual,residual,...
                  0,0,'Use the original Hamiltonian arithmetic, not a complex zero-copy.');
result=test_true(result,'Exact-reference validity is separate from budget',...
                   all(d.reduction.valid,'all')&&any(d.reduction.abs_error_bound>...
                   p.cce.numerics.atol+p.cce.numerics.rtol*abs(d.reduction.value),'all'),...
                   'A finite complete-bath density remains valid above the requested budget.');
result=test_true(result,'Lower-order masks retain propagated errors',...
                   any(~d.reduction.by_order_valid(:,:,1),'all'),...
                   'The lower-order numerical budget still masks unresolved elements.');

% Projected spectral products retain arbitrary preparation and receiver scaling
sys=struct('magnet',0.03,'isotopes',{{'E','1H'}},...
           'output','hush','parallel',{{'local',1}});
inter=struct(); inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2); inter.coupling.matrix{1,2}=[0 0 2;0 0 1;0 0 3];
system=create(sys,inter); p=struct('observed',1,'times',[0 1e-5 3e-5]);
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce=struct('approximation','projected','order',1);
p.preparation=struct('states',{{[1 -2i;2i 4]/5}},'spins',1);
p.coil=exp(0.37i)*[0 1;0 0]; p.normalized=false;
[~,diagnostics]=clexpn(system,@hahn_echo,p,'labframe');
model=cce_compile(system,diagnostics.parameters,'labframe');
config=rmfield(diagnostics.ensemble.spin_config,'assignment');
[cluster,~]=cce_subsystem(model,2,config,[]);
[H,R,K,p]=cce_generators(cluster,diagnostics.parameters,[]);
[actual,aux]=hahn_echo(cluster,p,H,R,K);
plain=rmfield(H,{'hamiltonian_abs_error_bound','hamiltonian_norm'});
expected=hahn_echo(cluster,p,plain,R,K);
result=test_close(result,'Projected unnormalised spectral pulse products',...
                  actual,expected,2e-14,2e-14,...
                  'Arbitrary coherent preparation and receiver phase retain their scale.');
result=test_true(result,'Projected spectral diagnostics',...
                 all(aux.abs_error_bound>0)&&aux.hamiltonian_bound>=0,...
                 'Numerical bounds and backward residuals remain attached to the signal.');
p.normalized=true;
actual=hahn_echo(cluster,p,H,R,K); expected=hahn_echo(cluster,p,plain,R,K);
result=test_close(result,'Projected normalized spectral pulse products',...
                  actual,expected,2e-14,2e-14,...
                  'Normalization uses the detected initial amplitude, including its phase.');

% A six-state bath retains the evaluated mixed-density trace at zero time
sys.magnet=0; sys.isotopes={'E','E','14N'};
inter=struct(); inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
system=create(sys,inter); p=struct();
p.observed=1; p.transition=struct('spins',1,'ms',[1/2 -1/2]);
p.times=0; p.cce.approximation='projected'; p.cce.order=1;
p.cce.expansion='partition'; p.cce.groups={[2 3]};
p.cce.partition=struct('mode','user','max_dimension',6);
p.cce.graph=struct('mode','complete');
actual=clexpn(system,@hahn_echo,p,'labframe');
rho_b=kron(eye(2)/2,eye(3)/3);
result=test_close(result,'Projected six-state zero-time trace',...
                  actual,trace(rho_b),0,0,...
                  'Zero-duration propagation evaluates the supplied density trace.');

end


