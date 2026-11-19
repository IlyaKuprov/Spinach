% Validate closed-system experiments against independent matrix propagation.
%
% Syntax:
%
%                      result=test_clexpn_experiments()
%
% Outputs:
%
%     result - standard registered-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_experiments()

% Noncommuting bath terms produce a complex, nontrivial echo modulation
result=new_test_result('kernel/clexpn_experiments','Native cluster experiments',...
                       'Free evolution, Hahn, CPMG, ESEEM and spectral bath pumping.');
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
inter.zeeman.matrix={diag([0 0 3]),diag([0 0 11])};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;5 0 7];
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.007 0.019];
p.transition.spins=1;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
p.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8);
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2])}});
[p,~]=cce_options(spin_system,p,@free_evolution,'labframe');
model=cce_compile(spin_system,p,'labframe');
[cluster,~]=cce_subsystem(model,2,struct(),[]);
[H,R,K,p]=cce_generators(cluster,p,[]);
[rho0,coil,reference]=cce_observation(cluster,p);
pulse=kron([0 -1i;-1i 0],eye(2));
sequences={@free_evolution,@hahn_echo,@cpmg_echo,@cce_eseem};
p.nloops=3;
fractions={[],0.5,([1 2 3]-0.5)/3,0.5};
for j=1:numel(sequences)

    % Compare each named sequence with dense propagation
    actual=sequences{j}(cluster,p,H,R,K);
    expected=zeros(size(actual));
    for n=1:numel(p.times)

        % Propagate and detect the current observation time
        u=eye(4);
        previous=0;
        if p.times(n)>0
            for fraction=fractions{j}

                % Apply the next free interval and physical rotation
                u=pulse*expm(-1i*H*(fraction-previous)*p.times(n))*u;
                previous=fraction;

            end
        end
        u=expm(-1i*H*(1-previous)*p.times(n))*u;
        expected(n)=trace(coil'*u*rho0*u')/reference;

    end
    result=test_close(result,['Independent ' func2str(sequences{j})],actual,expected,3e-13,3e-13,...
                      'Dense expm independently checks elapsed time, pulses and complex receiver.');

end

% Complex preparation exposes odd projected refocusing phases in both representations
phased=p;
phased.rho0=[0.5 -0.5i;0.5i 0.5];
phased.normalized=true;
phased.coil=exp(0.37i)*[0 1;0 0];
invariant=spin_system;
invariant.inter.coupling.matrix{1,2}=[0 0 5;0 0 0;0 0 7];
invariant.inter.coupling.matrix{2,1}=zeros(3);
invariant.inter.zeeman.matrix{2}=diag([0 0 11]);
phase_model=cce_compile(invariant,rmfield(phased,'cce_aux'),'labframe');
[phase_cluster,~]=cce_subsystem(phase_model,2,struct(),[]);
for sequence={@cpmg_echo,@cce_eseem}

    % Select the refocusing experiment
    for count=[1 2 3]

        % Propagate and detect the current observation time
        phased.nloops=count;
        for formalism={'zeeman-hilb','zeeman-liouv'}

            % Compare Hilbert and Liouville representations
            full_cluster=phase_cluster;
            full_cluster.formalism=formalism{1};
            phased.cce.approximation='generalized';
            [hf,rf,kf,full_p]=cce_generators(full_cluster,phased,[]);
            expected=sequence{1}(full_cluster,full_p,hf,rf,kf);
            phased.cce.approximation='projected';
            [hp,rp,kp,proj_p]=cce_generators(full_cluster,phased,[]);
            actual=sequence{1}(full_cluster,proj_p,hp,rp,kp);
            result=test_close(result,'Projected fixture manifold invariance',...
                              proj_p.cce_aux.projection_leakage,0,0,0,...
                              'The comparison fixture has no observed-spin population transfer.');
            result=test_close(result,sprintf('Projected phase %s %d %s',...
                              func2str(sequence{1}),count,formalism{1}),...
                              actual,expected,3e-13,3e-13,...
                              'Physical coherence, complex receiver, and zero-time behaviour match full-space evolution.');

        end

    end

end

% A fixed receiver detects the physical coherence in every projected Hahn path
for sites={[],2}

    % Accumulate the retained physical-spin contribution
    [retained,~]=cce_subsystem(phase_model,sites{1},struct(),[]);
    for sequence={@hahn_echo,@bath_pumped_hahn}

        % Select the fixed-receiver experiment
        for formalism={'zeeman-hilb','zeeman-liouv'}

            % Compare Hilbert and Liouville representations
            retained.formalism=formalism{1};
            for normalised=[false true]

                % Compare normalised and unnormalised detection
                q=phased;
                q.normalized=normalised;
                if isequal(sequence{1},@bath_pumped_hahn)
                    q.pump=struct('center',11,'halfwidth',0,'units','rad/s',...
                                  'fraction',0.25,'allow_empty',true);
                end
                q.cce.approximation='generalized';
                [hf,rf,kf,full_p]=cce_generators(retained,q,[]);
                expected=sequence{1}(retained,full_p,hf,rf,kf);
                q.cce.approximation='projected';
                [hp,rp,kp,proj_p]=cce_generators(retained,q,[]);
                actual=sequence{1}(retained,proj_p,hp,rp,kp);
                result=test_close(result,'Hahn physical fixed receiver',actual,expected,...
                                  3e-13,3e-13,'Complex preparation and detector match full-space evolution.');
                if isequal(sequence{1},@hahn_echo)&&strcmp(formalism{1},'zeeman-hilb')
                    plain=rmfield(hp,{'hamiltonian_abs_error_bound','hamiltonian_norm'});
                    fallback=hahn_echo(retained,proj_p,plain,rp,kp);
                    result=test_close(result,'Hahn non-spectral fixed receiver',fallback,expected,...
                                      3e-13,3e-13,'The interval fallback detects the same physical coherence.');
                end

            end

        end

    end

end

% Nonzero dissipation and manifold-preserving rotations retain the fixed receiver
[retained,~]=cce_subsystem(phase_model,2,struct(),[]);
retained.formalism='zeeman-liouv';
q=phased;
q.bath_jumps=struct('support',2,'operator','L+','rate',2);
for angle=[pi 2*pi]

    % Set an interchanging or preserving rotation
    q.refocusing=struct('spins',1,'operator',angle*[0 1;1 0]/2,'time',0.007);
    q.cce.approximation='generalized';
    [hf,rf,kf,full_p]=cce_generators(retained,q,[]);
    expected=hahn_echo(retained,full_p,hf,rf,kf);
    q.cce.approximation='projected';
    [hp,rp,kp,proj_p]=cce_generators(retained,q,[]);
    actual=hahn_echo(retained,proj_p,hp,rp,kp);
    result=test_close(result,'Dissipative Hahn fixed receiver',actual,expected,...
                      3e-13,3e-13,'Applied and absent swaps use the original receiver and reference.');

end

% Default preparation also checks the bath-pumped spectral path against full space
q=rmfield(phased,{'rho0','coil'});
q.normalized=true;
q.pump=struct('center',11,'halfwidth',0,'units','rad/s',...
              'fraction',0.25,'allow_empty',false);
retained.formalism='zeeman-hilb';
q.cce.approximation='generalized';
[hf,rf,kf,full_p]=cce_generators(retained,q,[]);
expected=bath_pumped_hahn(retained,full_p,hf,rf,kf);
q.cce.approximation='projected';
[hp,rp,kp,proj_p]=cce_generators(retained,q,[]);
actual=bath_pumped_hahn(retained,proj_p,hp,rp,kp);
result=test_close(result,'Bath-pumped spectral fixed receiver',actual,expected,...
                  3e-13,3e-13,'Spectral propagation retains the full-space physical receiver.');

% Empty spectral selection is a documented no-pump control
p.pump=struct('center',1e12,'halfwidth',0,'units','Hz','fraction',0.25,'allow_empty',true);
actual=bath_pumped_hahn(cluster,p,H,R,K);
expected=hahn_echo(cluster,p,H,R,K);
result=test_close(result,'Empty allowed pump window',actual,expected,3e-13,3e-13,...
                  'An empty explicit window leaves the Hahn experiment unchanged.');

% Free-evolution density has the same column-major ordering as direct partial trace
p.output='density';
actual=free_evolution(cluster,p,H,R,K);
expected=zeros(3,4);
for n=1:numel(p.times)

    % Propagate and detect the current observation time
    u=expm(-1i*H*p.times(n));
    rho=u*rho0*u';
    observed=zeros(2);
    for a=1:2

        % Trace each observed-spin row over the bath
        for b=1:2
            observed(a,b)=trace(rho((a-1)*2+(1:2),(b-1)*2+(1:2)));
        end

    end
    expected(n,:)=observed(:).';

end
result=test_close(result,'Column-major reduced density',actual,expected,3e-13,3e-13,...
                  'Direct free evolution returns physical observed density elements.');

% A physical receiver phase remains complex before normalisation
phase=0.37;
phased=p;
phased.output='observable';
phased.normalized=false;
phased.coil=exp(1i*phase)*model.transition.alpha*model.transition.beta';
unphased=phased;
unphased.coil=model.transition.alpha*model.transition.beta';
actual=free_evolution(cluster,phased,H,R,K);
expected=exp(-1i*phase)*free_evolution(cluster,unphased,H,R,K);
result=test_close(result,'Complex receiver phase',actual,expected,3e-13,3e-13,...
                  'Detection uses the conjugate-transposed receiver operator.');

% Inclusive spectral-window boundaries select the same physical bath pulse
p.output='signal';
p.pump.center=model.zeeman{2}(3,3);
p.pump.units='rad/s';
p.pump.halfwidth=0;
p.pump.allow_empty=false;
actual=bath_pumped_hahn(cluster,p,H,R,K);
reference_pump=p;
reference_pump.pump.halfwidth=eps(abs(p.pump.center));
expected=bath_pumped_hahn(cluster,reference_pump,H,R,K);
result=test_close(result,'Inclusive pump boundary',actual,expected,0,0,...
                  'A site exactly on the spectral boundary is retained.');

% Finite drives remain experiment-owned and retain the coupled drift
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
system=create(sys,inter);
p=struct();
p.observed=1;
p.times=[0 0.2];
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.output='density';
p.cce=struct('approximation','generalized','order',1,'graph',struct('mode','complete'));
p.drive=struct('spins',1,'axis','z','time',0,...
               'duration',0.2,'amplitude',7,'detuning',0);
actual=clexpn(system,@finite_drive,p,'labframe');
rho=[1 1;1 1]/2;
U=expm(-1i*7*diag([0.5 -0.5])*0.2);
expected=U*rho*U';
result=test_close(result,'Finite-drive acceptance probe',actual(2,:),expected(:).',...
                  2e-14,0,'The native finite drive rotates the observed density.');
p.preparation='transverse';
named=clexpn(system,@finite_drive,p,'labframe');
result=test_close(result,'Transverse preparation',named,actual,2e-14,0,...
                  'The named positive preparation matches the transition superposition.');
p.cce.basis.formalism='zeeman-liouv';
liouville=clexpn(system,@finite_drive,p,'labframe');
result=test_close(result,'Finite-drive Liouville representation',liouville,actual,...
                  2e-14,0,'The Hamiltonian commutator gives the same observed density.');
p.cce=rmfield(p.cce,'basis');
p.times=[0 0.05 0.2 0.4];
p.drive.duration=0.1;
p.drive.time=0.05;
p.drive.detuning=2;
p.drive(2)=p.drive(1);
p.drive(2).time=0.1;
p.drive(2).axis='x';
p.drive(2).amplitude=3;
actual=clexpn(system,@finite_drive,p,'labframe');
expected=zeros(size(actual));
Sz=diag([0.5 -0.5]);
Sx=[0 1;1 0]/2;
for n=1:numel(p.times)

    % Propagate and detect the current observation time
    t=p.times(n);
    state=rho;
    cuts=unique([0 t min([0.05 0.1 0.15 0.2],t)]);
    for k=1:numel(cuts)-1

        % Propagate the simultaneous active drive intervals
        middle=(cuts(k)+cuts(k+1))/2;
        generator=zeros(2);
        if middle>0.05&&middle<0.15
            generator=generator+9*Sz;
        end
        if middle>0.1&&middle<0.2
            generator=generator+3*Sx+2*Sz;
        end
        U=expm(-1i*generator*(cuts(k+1)-cuts(k)));
        state=U*state*U';

    end
    expected(n,:)=state(:).';

end
result=test_close(result,'Overlapping finite intervals',actual,expected,3e-14,0,...
                  'Start/stop boundaries, detuning and partial intervals match dense expm.');
p.drive=p.drive(1);
p.drive.axis=[0 0 1];
p.drive=rmfield(p.drive,'time');
p.drive.fraction=0.5;
actual=clexpn(system,@finite_drive,p,'labframe');
for n=1:numel(p.times)

    % Propagate and detect the current observation time
    U=expm(-1i*9*Sz*min(0.1,p.times(n)/2));
    state=U*rho*U';
    expected(n,:)=state(:).';

end
result=test_close(result,'Fractional finite-drive placement',actual,expected,3e-14,0,...
                  'The duration remains in seconds while start scales with elapsed time.');

% Noncommuting drift remains active throughout a finite control interval
inter.coupling.matrix{1,2}=diag([1 2 3]);
system=create(sys,inter);
p.drive=rmfield(p.drive,'fraction');
p.drive.time=0;
p.drive.duration=0.1;
p.drive.axis='x';
p.preparation='transverse';
actual=clexpn(system,@finite_drive,p,'labframe');
Sy=[0 -1i;1i 0]/2;
drift=2*pi*(kron(Sx,Sx)+2*kron(Sy,Sy)+3*kron(Sz,Sz));
drive=kron(7*Sx+2*Sz,eye(2));
for n=1:numel(p.times)

    % Propagate and detect the current observation time
    duration=min(0.1,p.times(n));
    U=expm(-1i*drift*(p.times(n)-duration))*expm(-1i*(drift+drive)*duration);
    state=U*kron(rho,eye(2)/2)*U';
    observed=zeros(2);
    for a=1:2

        % Trace each observed-spin row over the bath
        for b=1:2
            observed(a,b)=trace(state((a-1)*2+(1:2),(b-1)*2+(1:2)));
        end

    end
    expected(n,:)=observed(:).';

end
result=test_close(result,'Finite-drive coupled drift',actual,expected,6e-14,0,...
                  'Dense full-bath propagation includes noncommuting dipolar drift.');

% Named finite drives compile the same control without a public event array
q=rmfield(p,'drive');
q.drive=struct('spins',1,'axis','x','amplitude',7,'detuning',2,'duration',0.1);
selected=p.times==0|p.times>=q.drive.duration;
q.times=p.times(selected);
named=clexpn(system,@finite_drive,q,'labframe');
result=test_close(result,'Named finite-drive contract',named,actual(selected,:),0,0,...
                  'The declarative drive preserves the controlled density bitwise.');
invalid={struct('spins',2),struct('axis',[1 1 0]),struct('amplitude',Inf),...
         struct('detuning',1i),struct('duration',0),struct('extra',1)};
for n=1:numel(invalid)

    % Reject the invalid drive descriptor
    bad=q;
    field=fieldnames(invalid{n});
    bad.drive.(field{1})=invalid{n}.(field{1});
    rejected=false;
    try
        clexpn(system,@finite_drive,bad,'labframe');
    catch
        rejected=true;
    end
    result=test_true(result,['Named finite-drive validation ' field{1}],rejected,...
                      'Every declarative drive field is validated.');

end
bad=q;
bad.times=[0 0.02];
partial=clexpn(system,@finite_drive,bad,'labframe');
U=expm(-1i*(drift+drive)*0.02);
state=U*kron(rho,eye(2)/2)*U';
observed=zeros(2);
for a=1:2

    % Trace each observed-spin row over the bath
    for b=1:2
        observed(a,b)=trace(state((a-1)*2+(1:2),(b-1)*2+(1:2)));
    end

end
result=test_close(result,'Named finite-drive observation time',partial(2,:),observed(:).',...
                  6e-14,0,'An observation during a drive retains simultaneous coupled evolution.');
bad=q;
bad.sequence=struct('type','finite','target',1,'axis','x',...
                    'placement','absolute','at',0,'duration',0.1,'amplitude',7,'detuning',2);
rejected=false;
try
    clexpn(system,@finite_drive,bad,'labframe');
catch
    rejected=true;
end
result=test_true(result,'Named finite-drive event ownership',rejected,...
                  'A competing event array is rejected.');

% Mixed finite and ideal controls retain their explicit representation limit
bad=q;
bad.cce.basis.formalism='zeeman-liouv';
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
inter.coupling.matrix=cell(2);
system=create(sys,inter);
p.cce.basis.formalism='zeeman-liouv';
p.bath_jumps=struct('support',1,'operator','Lz','rate',2);
actual=clexpn(system,@finite_drive,p,'labframe');
J=sqrt(2)*Sz;
loss=J'*J;
D=kron(conj(J),J)-(kron(eye(2),loss)+kron(loss.',eye(2)))/2;
A=-1i*(kron(eye(2),7*Sx+2*Sz)-kron((7*Sx+2*Sz).',eye(2)))+D;
for n=1:numel(p.times)

    % Propagate and detect the current observation time
    duration=min(0.1,p.times(n));
    expected(n,:)=(expm(D*(p.times(n)-duration))*expm(A*duration)*rho(:)).';

end
result=test_close(result,'Finite-drive GKSL evolution',actual,expected,6e-14,0,...
                  'The dissipator remains active during and after the drive.');
p=rmfield(p,'bath_jumps');
p.cce=rmfield(p.cce,'basis');

% Raw finite, ideal, and empty event arrays are all removed inputs
for events={struct('type','finite'),struct('type','ideal'),struct([]),[]}

    % Reject the removed raw sequence representation
    rejected=p;
    rejected=rmfield(rejected,'drive');
    rejected.sequence=events{1};
    failed=false;
    try
        clexpn(system,@free_evolution,rejected,'labframe');
    catch err
        failed=strcmp(err.identifier,'cce_options:sequence');

    end
    result=test_true(result,'Free evolution rejects raw events',failed,...
                     'Use finite_drive, refocusing_train, or a named echo experiment.');

end
p.cce.approximation='projected';
p.output='signal';
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
positions={[2 3],[1 3],[1 4],[1 9]};
dimensions=[4 3 4 16];
for k=1:numel(isotopes)

    % Prepare the specified product spin manifolds
    sys.isotopes=isotopes{k};
    inter=struct();
    inter.zeeman.matrix=repmat({zeros(3)},1,numel(sys.isotopes));
    inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
    system=create(sys,inter);
    p=struct();
    p.observed=1:numel(sys.isotopes);
    p.times=0;
    p.output='density';
    p.transition=struct('spins',p.observed,'ms',selectors{k});
    p.cce=struct('approximation','generalized','order',0);
    actual=clexpn(system,@free_evolution,p,'labframe');
    vector=zeros(dimensions(k),1);
    vector(positions{k})=1/sqrt(2);
    expected=vector*vector';
    result=test_close(result,sprintf('Example preparation selector %d',k),...
                      actual,expected(:).',2e-14,0,...
                      'Selected product manifolds replace hand-built preparation vectors.');

end

% Coherent diagnostics retain the original arithmetic and exact-reference masks
sys.magnet=0;
sys.isotopes={'E','1H'};
sys.output='hush';
inter=struct();
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([1 2 3]);
system=create(sys,inter);
p=struct();
p.observed=1;
p.transition.spins=1;
p.transition.ms=[0.5 -0.5];
p.times=[0 0.1];
p.output='density';
p.cce.approximation='generalized';
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
configuration=d.ensemble.spin_config;
configuration=rmfield(configuration,'assignment');
[retained,~]=cce_subsystem(model,2,configuration,[]);
[H,R,K,pseq]=cce_generators(retained,d.parameters,[]);
[~,aux]=free_evolution(retained,pseq,H,R,K);
[v,e]=eig(full((H+H')/2),'vector');
residual=norm(H*v-v.*e.','fro');
orthogonality=norm(v'*v-eye(size(H)),'fro');
operations=12*size(H,1)^3;
gamma=operations*eps/(1-operations*eps);
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
inter=struct();
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=[0 0 2;0 0 1;0 0 3];
system=create(sys,inter);
p=struct('observed',1,'times',[0 1e-5 3e-5]);
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce=struct('approximation','projected','order',1);
p.preparation=struct('states',{{[1 -2i;2i 4]/5}},'spins',1);
p.coil=exp(0.37i)*[0 1;0 0];
p.normalized=false;
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
actual=hahn_echo(cluster,p,H,R,K);
expected=hahn_echo(cluster,p,plain,R,K);
result=test_close(result,'Projected normalized spectral pulse products',...
                  actual,expected,2e-14,2e-14,...
                  'Normalization uses the detected initial amplitude, including its phase.');

% A six-state bath retains the evaluated mixed-density trace at zero time
sys.magnet=0;
sys.isotopes={'E','E','14N'};
inter=struct();
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
system=create(sys,inter);
p=struct();
p.observed=1;
p.transition=struct('spins',1,'ms',[1/2 -1/2]);
p.times=0;
p.cce.approximation='projected';
p.cce.order=1;
p.cce.expansion='partition';
p.cce.groups={[2 3]};
p.cce.partition=struct('mode','user','max_dimension',6);
p.cce.graph=struct('mode','complete');
actual=clexpn(system,@hahn_echo,p,'labframe');
rho_b=kron(eye(2)/2,eye(3)/3);
result=test_close(result,'Projected six-state zero-time trace',...
                  actual,trace(rho_b),0,0,...
                  'Zero-duration propagation evaluates the supplied density trace.');

% A complete bath can yield a valid echo above a strict propagation budget
q=struct('observed',1,'times',[0 0.019],'output','observable','nloops',3);
q.transition=struct('spins',1,'ms',[0.5 -0.5]);
q.rho0=[0.5 -0.4i;0.4i 0.5];
detector=3*exp(0.37i)*[0 1;0 0];
q.cce.order=1;
q.cce.graph.mode='complete';
q.cce.numerics=struct('atol',1e-30,'rtol',1e-30,'max_memory_bytes',1e8);
for sequence={@cpmg_echo,@cce_eseem,@hahn_echo,@bath_pumped_hahn}

    % Count the nonzero free intervals separated by ideal inversions
    intervals=2;
    if isequal(sequence{1},@cpmg_echo)
        intervals=q.nloops+1;
    end
    hahn=isequal(sequence{1},@hahn_echo);
    pumped=isequal(sequence{1},@bath_pumped_hahn);
    formalisms={'zeeman-hilb','zeeman-liouv'};
    if hahn
        formalisms={'zeeman-hilb'};
    end
    if hahn||pumped
        q.pump=struct('center',11,'halfwidth',0,'units','rad/s','fraction',0.5,'allow_empty',false);
    end
    for approximation={'projected','generalized'}

        % Generalised detection supports independently scaled complex receivers
        q.cce.approximation=approximation{1};
        q.coil=detector;
        if strcmp(approximation{1},'generalized')
            q.coil={detector;2i*detector};
        end
        for formalism=formalisms

            % Compare the segment sum with the norm-preserving closed-system estimate
            q.cce.basis.formalism=formalism{1};
            for normalised=[false true]

                % Exercise the production dispatch, exact reduction, and budget decision
                q.normalized=normalised;
                [~,d]=clexpn(invariant,sequence{1},q,'labframe');
                model=cce_compile(invariant,d.parameters,'labframe');
                model.formalism=formalism{1};
                checks=cce_diagnostics(model,d.raw,d.reduction,d.parameters);
                [retained,~]=cce_subsystem(model,2,struct(),[]);
                [H,R,K,local]=cce_generators(retained,d.parameters,[]);
                [rho,coils,reference]=cce_observation(retained,local);
                if ~iscell(coils)
                    coils={coils};
                end
                if isstruct(H)&&~isempty(H.dissipator)
                    db=retained.bath_dimension;
                    scale=64*eps*db^2;
                    frequency=norm(-1i*(kron(eye(db),H.left)-kron(H.right.',eye(db)))+H.dissipator,1);
                elseif isstruct(H)
                    scale=16*eps*retained.bath_dimension^3;
                    frequency=norm(H.left,1)+norm(H.right,1);
                elseif strcmp(formalism{1},'zeeman-hilb')
                    scale=16*eps*retained.dimension^3;
                    frequency=2*norm(H,1);
                else
                    scale=64*eps*size(H,1);
                    frequency=norm(H+1i*R+1i*K,1);
                end
                weights=cellfun(@(x)norm(x,'fro'),coils(:).');
                if normalised
                    weights=weights./abs(reference);
                end
                pulse_cost=0;
                if hahn||pumped
                    dim=retained.dimension;
                    products=1+pumped;
                    if isstruct(H)
                        dim=retained.bath_dimension;
                        products=pumped;
                    end
                    pulse_cost=products*16*dim^3*eps;
                end
                expected=(scale*(intervals+q.times(2)*frequency)+pulse_cost)*norm(rho,'fro')*weights;
                label=sprintf('%s %s %s normalised=%d',func2str(sequence{1}),...
                              approximation{1},formalism{1},normalised);
                result=test_close(result,['Echo segment estimates ' label],...
                                  d.reduction.abs_error_bound(2,:)./expected,ones(size(expected)),0,5e-14,...
                                  'Ideal pulses preserve the state norm; detection scales each segment sum.');
                result=test_true(result,['Echo strict budget ' label],...
                                 all(checks.valid,'all')&&all(d.reduction.abs_error_bound(1,:)==0)&&...
                                 all(checks.numerical_budget_met(1,:))&&~any(checks.numerical_budget_met(2,:)),...
                                 'Only the unevolved signal meets the strict numerical budget.');

            end

        end

    end

end

% A relaxing bath can amplify Frobenius error without amplifying trace norm
sys.isotopes={'1H','13C'};
system=create(sys,struct());
q.cce.basis.formalism='zeeman-liouv';
q.pump.center=0;
q.times=[0 1];
q.bath_jumps=struct('support',2,'operator','L+','rate',4);
q.normalized=false;
q.coil=detector;
q.output='signal';
for approximation={'projected','generalized'}

    % Compare local errors against the analytic polarisation of an initially mixed bath
    q.cce.approximation=approximation{1};
    model=cce_compile(system,q,'labframe');
    model.formalism='zeeman-liouv';
    [retained,~]=cce_subsystem(model,2,struct(),[]);
    [H,R,K,local]=cce_generators(retained,model.parameters,[]);
    [~,aux]=bath_pumped_hahn(retained,local,H,R,K);
    if strcmp(approximation{1},'projected')
        dim=2;
        initial=abs(q.rho0(1,2));
        products=1;
    else
        dim=4;
        initial=norm(q.rho0,'fro');
        products=2;
    end
    half=exp(-q.bath_jumps.rate/2)/2;
    before=initial/sqrt(2);
    during=initial*norm([1-half half]);
    propagation=64*eps*dim^2*(1+q.bath_jumps.rate);
    expected=sqrt(dim)*sqrt(2)*norm(detector,'fro')*...
             (propagation*(before+during)+products*16*dim^3*eps*during);
    result=test_close(result,['Relaxing bath-pump estimate ' approximation{1}],...
                      aux.abs_error_bound(2)/expected,1,0,5e-14,...
                      'Both half intervals, simultaneous pulse products, and the channel gain are retained.');

end

% Each partial-trace component has the same Frobenius sensitivity
q=rmfield(q,{'bath_jumps','coil'});
q.cce.basis.formalism='zeeman-hilb';
q.output='density';
[~,d]=clexpn(system,@hahn_echo,q,'labframe');
expected=sqrt(2)*3*16*4^3*eps*norm(kron(q.rho0,eye(2)/2),'fro');
result=test_close(result,'Interval Hahn partial-trace estimate',...
                  d.raw{1}(end).abs_error_bound(2,:)/expected,ones(1,4),0,5e-14,...
                  'The two free intervals and the observed inversion contribute before partial tracing.');

end


