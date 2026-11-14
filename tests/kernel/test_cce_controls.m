% Finite controls, named timestamps and exterior-field regression tests.
% Syntax: result=test_cce_controls()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_controls()

% Initialise the standard harness result and a driven two-spin system
result=new_test_result('kernel/cce_controls','CCE realistic controls',...
                       'Check finite drives, waveform refinement and named pulse schedules.');
[model,p,ensemble]=fixture(); spin=pauli(2); x=spin.x; z=spin.z; context=ensemble.contexts(1);
p.sequence=struct('type','finite','target',1,'placement','absolute','at',0,...
                  'duration',0.1,'axis','x','amplitude',7,'detuning',3);
s=cce_controls(model,2,context,[],p); H=cce_hamiltonian(model,2,struct());
U=expm(-1i*0.1*(H+7*kron(x,eye(2))+3*kron(z,eye(2))));
initial=kron(context.central_rho,context.bath_rho{1}); rho=U*initial*U'; expected=partial(rho);
error=norm(s.raw.density(:,:,end)-expected,'fro');
result=test_close(result,'Finite detuned drive matches independent dense expm',s.raw.density(:,:,end),expected,1e-12,0,'Finite detuned drive matches independent dense expm');

% Require no-op control dispatch to preserve explicit and default normalization
q=p; q.quantity='coherence'; q.central.alpha=[1;0]; q.central.beta=[0;1];
q.central.rho=ones(2)/2; q.sequence=[]; ss=model.spin_system;
for normalized=[false true]
    q.normalized=normalized; plain=cce(ss,q);
    q.sequence=struct('type','finite','target',1,'placement','absolute','at',0,...
                      'duration',0.01,'axis','x','amplitude',0,'detuning',0);
    driven=cce(ss,q);
    result=test_close(result,'Control dispatch preserves coherence normalization',...
                      driven.value,plain.value,1e-12,1e-12,...
                      'A zero-amplitude control must not change the requested coherence scale.');
    q.sequence=[];
end
q=rmfield(q,'normalized'); plain=cce(ss,q);
q.sequence=struct('type','finite','target',1,'placement','absolute','at',0,...
                  'duration',0.01,'axis','x','amplitude',0,'detuning',0);
driven=cce(ss,q);
result=test_close(result,'Control dispatch preserves default raw coherence',...
                  driven.value,plain.value,1e-12,1e-12,...
                  'Omitting normalized must retain the same raw coherence on both dispatch paths.');
fprintf('Control normalization dispatch discrepancy %.3g.\n',max(abs(driven.value-plain.value)));

% Compare named and explicit Hahn dispatch with the same physical events
q.sequence=struct('type','ideal','target',1,'placement','fraction',...
                  'at',0.5,'axis','x','angle',pi);
explicit_echo=cce(ss,q);
q.sequence=struct('name','hahn','count',1,'axis','x','target',1);
named_echo=cce(ss,q);
result=test_close(result,'Named and explicit Hahn normalization',...
                  named_echo.value,explicit_echo.value,1e-12,1e-12,...
                  'Naming the same pulse sequence must not change the default coherence scale.');

% Exercise the public driver through the named/finite dispatch hooks
integrated=cce(model.spin_system,p);
result=test_close(result,'Public finite-control driver',reshape(integrated.value(end,:),2,2),expected,...
                  1e-12,0,'Driver consumes the identical segment propagator');

% Approach the ideal pulse with a fixed rotation angle and decreasing width
ideal=struct('type','ideal','target',1,'placement','absolute','at',0.05,'axis','x','angle',pi);
p.sequence=ideal; reference=cce_controls(model,2,context,[],p); errors=zeros(1,3);
for k=1:3
    width=10^(-k-2); p.sequence=struct('type','finite','target',1,'placement','absolute',...
        'at',0.05-width/2,'duration',width,'axis','x','amplitude',pi/width,'detuning',0);
    s=cce_controls(model,2,context,[],p); errors(k)=norm(s.raw.density(:,:,end)-reference.raw.density(:,:,end),'fro');
end
result=test_true(result,'Finite pulse approaches ideal rotation',all(diff(errors)<0)&&errors(end)<2e-5,'Finite pulse approaches ideal rotation');

% Compare midpoint piecewise envelopes against a much finer independent product
mesh=[8 16 32]; waveform_errors=zeros(size(mesh)); reference=waveform(model,context,1024);
for k=1:numel(mesh)
    n=mesh(k); dt=0.1/n;
    p.sequence=repmat(struct('type','finite','target',1,'placement','absolute','at',0,...
                    'duration',dt,'axis','x','amplitude',0,'detuning',0),1,n);
    for j=1:n, p.sequence(j).at=(j-1)*dt; p.sequence(j).amplitude=10*sin(pi*(j-0.5)/n); end
    s=cce_controls(model,2,context,[],p); waveform_errors(k)=norm(s.raw.density(:,:,end)-reference,'fro');
end
result=test_true(result,'Piecewise midpoint refinement is second order',all(waveform_errors(1:end-1)./waveform_errors(2:end)>3.8),'Piecewise midpoint refinement is second order');

% Check literal XY8 phases and analytic UDD timestamps
p.sequence=struct('name','xy8','count',8,'axis','x','target',1);
s=cce_controls(model,2,[],[],p);
result=test_close(result,'XY8 timestamps',s.schedule.positions(:,end),0.1*((1:8)'-0.5)/8,1e-15,0,'XY8 timestamps');
result=test_true(result,'XY8 phases',strcmp([s.schedule.events.axis],'xyxyyxyx'),'XY8 phases');
p.sequence=struct('name','udd','count',5,'axis','y','target',1); s=cce_controls(model,2,[],[],p);
result=test_close(result,'UDD timestamps',s.schedule.positions(:,end),0.1*sin(pi*(1:5)'/12).^2,1e-15,0,'UDD timestamps');

% Exercise time-dependent exterior mean fields and explicit finite-drive refusal
p.meanfield='static'; p.sequence=ideal; p.sequence.target=2;
s=cce_controls(model,[],context,[],p); field=s.exterior_expectations{end};
result=test_close(result,'Exterior ideal pulse flips the sampled field',field{2}(1,3),-context.expectations(1,3),1e-14,0,'Exterior ideal pulse flips the sampled field');
Hbefore=s.hamiltonians{end}{2}; Hafter=s.hamiltonians{end}{3};
result=test_true(result,'Exterior pulse changes subsequent Hamiltonian',norm(Hbefore-Hafter,'fro')>1,'Exterior pulse changes subsequent Hamiltonian');
p.sequence=struct('type','finite','target',2,'placement','absolute','at',0,...
                  'duration',0.1,'axis','x','amplitude',7,'detuning',0);
refused=false;
try
    cce_controls(model,[],context,[],p); catch exception, refused=strcmp(exception.identifier,'cce_controls:exterior_finite'); end
result=test_true(result,'Finite exterior drive is refused rather than frozen',refused,'Finite exterior drive is refused rather than frozen');

% Refuse finite central controls in the conditional-branch approximation
p.meanfield='none'; p.method='cce'; p.sequence.target=1; refused=false;
try
    cce_controls(model,2,context,[],p); catch exception, refused=strcmp(exception.identifier,'cce_controls:conditional_drive'); end
result=test_true(result,'Conditional finite central drive refused',refused,'General central drives need gCCE');

% Verify an explicitly ordered register pulse against a dense two-spin operator
p.method='gcce'; p.sequence=struct('type','ideal','target',[2 1],'placement','absolute',...
    'at',0,'operator',kron(z,x),'angle',0.7);
s=cce_controls(model,2,context,[],p); expected_pulse=expm(-1i*0.7*kron(x,z));
result=test_close(result,'Selective register operator ordering',s.pulses{end}{1},expected_pulse,...
                  1e-13,0,'The target tensor order is remapped to canonical physical order');

% Distinguish a physical spin-one rotation from a selective two-level pulse
sys.isotopes={'14N'}; sys.magnet=0; sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.output='hush';
inter=struct(); spin_system=create(sys,inter); q=p; q.central_spins=1; q.bath_spins=[];
q.central.rho=diag([1 0 0]); q.order=0; q.bath_state=struct('mode','mixed_product');
q.sequence=struct('type','ideal','target',1,'placement','absolute','at',0,'axis','x','angle',pi);
spin_one=cce_prepare(spin_system,q); states=cce_states(spin_one,q);
physical=cce_controls(spin_one,[],states.contexts(1),[],q);
q.sequence.operator=[0 1 0;1 0 0;0 0 0]/2;
selective=cce_controls(spin_one,[],states.contexts(1),[],q);
result=test_close(result,'Physical spin-one pi pulse',physical.raw.density(:,:,end),diag([0 0 1]),...
                  1e-12,0,'A full spin-one rotation maps m=+1 to m=-1');
result=test_close(result,'Selective two-level pi pulse',selective.raw.density(:,:,end),diag([0 1 0]),...
                  1e-12,0,'The declared two-level pulse maps m=+1 to m=0');
result=test_true(result,'Physical and selective pulses differ',...
    norm(physical.raw.density(:,:,end)-selective.raw.density(:,:,end),'fro')>1,...
    'A spin-one Cartesian rotation must not silently become a two-level pulse');

% Apply the declared angle to explicit generators on the static path
sys=struct(); sys.magnet=0; sys.isotopes={'E'}; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'processes',1};
inter=struct(); inter.zeeman.matrix={zeros(3)}; spin_system=create(sys,inter);
stat.central_spins=1; stat.bath_spins=[]; stat.method='gcce'; stat.expansion='standard';
stat.quantity='density'; stat.order=0; stat.central.rho=diag([1 0]); stat.assumptions='labframe';
stat.orientation=[0 0 0]; stat.interaction_cutoff=Inf; stat.graph=struct('mode','complete');
stat.meanfield='none'; stat.times=[0 1]; stat.bath_state=struct('mode','mixed_product');
stat.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
stat.sequence=struct('type','ideal','target',1,'placement','fraction','at',0.5,...
                     'operator',[0 1;1 0]/2,'angle',pi);
rotated=cce(spin_system,stat); stat.sequence.angle=0; unrotated=cce(spin_system,stat);
stat.sequence=rmfield(stat.sequence,'angle'); stat.sequence.operator=pi*[0 1;1 0]/2;
folded=cce(spin_system,stat);
result=test_close(result,'Static explicit generator with angle pi',rotated.density(:,:,end),diag([0 1]),...
                  1e-12,0,'An explicit Sx generator with angle pi inverts the central population');
result=test_close(result,'Static explicit generator with zero angle',unrotated.density(:,:,end),diag([1 0]),...
                  1e-12,0,'A zero angle must leave the central state unchanged');
result=test_close(result,'Static angle-folded generator',folded.density(:,:,end),rotated.density(:,:,end),...
                  1e-12,0,'A missing angle defaults to one, so pi*Sx equals Sx with angle pi');

% Print numerical evidence without adding nonstandard harness fields
fprintf('controls dense error %.3g; ideal-limit errors %s; waveform errors %s\n',error,mat2str(errors,6),mat2str(waveform_errors,6));

end

% Trace the independent bath index in the dense Kronecker reference
function rho=partial(full)
rho=zeros(2);
for a=1:2
    for b=1:2, rho(a,b)=trace(full((a-1)*2+(1:2),(b-1)*2+(1:2))); end
end
end

% Independently propagate a sampled sinusoidal drive with a fine mesh
function rho=waveform(model,context,n)
spin=pauli(2); x=spin.x; H=cce_hamiltonian(model,2,struct()); U=eye(4); dt=0.1/n;
for j=1:n, U=expm(-1i*dt*(H+10*sin(pi*(j-0.5)/n)*kron(x,eye(2))))*U; end
rho=partial(U*kron(context.central_rho,context.bath_rho{1})*U');
end

% Build a small noncommuting central-bath model in angular-frequency units
function [model,p,ensemble]=fixture()
sys.isotopes={'1H','1H'}; sys.magnet=0; sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.output='hush';
inter.coupling.matrix=cell(2); inter.coupling.matrix{1,2}=diag([2 3 5])/(2*pi);
spin_system=create(sys,inter); p.central_spins=1; p.bath_spins=2;
p.method='gcce'; p.expansion='standard'; p.quantity='density'; p.times=[0 0.1]; p.order=1;
p.central.rho=[1 1;1 1]/2; p.assumptions='labframe'; p.orientation=[0 0 0];
p.interaction_cutoff=Inf; p.graph=struct('mode','complete','score','coupling_norm');
p.bath_state=struct('mode','explicit_product','rho',{{diag([1 0])}});
p.meanfield='none'; p.sequence=[]; p.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
model=cce_prepare(spin_system,p); ensemble=cce_states(model,p);
end


