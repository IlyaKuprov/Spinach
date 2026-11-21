% Finite controls, named timestamps, and exterior-field regression tests. Syntax:
%
%                    result=test_cce_controls()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_controls()

% Initialise the standard harness result and a driven two-spin system
result=new_test_result('kernel/cce_controls','CCE realistic controls',...
                       'Check finite drives, waveform refinement and named pulse schedules.');
[model,p,ensemble]=fixture();
spin=pauli(2);
x=spin.x;
z=spin.z;
spin_config=ensemble.spin_configs(1);
p.drive=struct('spins',1,'duration',0.1,'axis','x','amplitude',7,'detuning',3);
[~,s_aux]=controlled(model,2,spin_config,p,@finite_drive);
local=cce_subsystem(model,2,spin_config,[]);
H=cce_generators(local,p,[]);
U=expm(-1i*0.1*(H+7*kron(x,eye(2))+3*kron(z,eye(2))));
initial=kron(p.rho0,spin_config.bath_rho{1});
rho=U*initial*U';
expected=partial(rho);
error=norm(s_aux.density(:,:,end)-expected,'fro');
result=test_close(result,'Finite detuned drive matches independent dense expm',s_aux.density(:,:,end),expected,1e-12,0,'Finite detuned drive matches independent dense expm');

% Normalise driven observables against the supplied observed spin configuration
q=p;
q.output='observable';
q.coil=diag([1 2]);
q.normalized=true;
supplied_config=ensemble.spin_configs(1);
supplied_config.rho0=diag([0.25 0.75]);
supplied_config.bath_rho={diag([0 1])};
supplied_config.expectations=[0 0 -0.5];
q.rho0=supplied_config.rho0;
supplied_config=rmfield(supplied_config,'rho0');
driven=controlled(model,2,supplied_config,q,@finite_drive);
h_observed=H+7*kron(x,eye(2))+3*kron(z,eye(2));
U=expm(-1i*h_observed*0.1);
rho=partial(U*kron(q.rho0,supplied_config.bath_rho{1})*U');
config_expected=[1;trace(q.coil*rho)/trace(q.coil*q.rho0)];
result=test_close(result,'Finite controls supplied spin_config reference',driven,config_expected,...
                  2e-12,2e-12,'Both the driven state and normalization use the supplied spin_config, not parameter rho.');
fprintf('R13 driven spin_config error %.17g initial %.17g\n',norm(driven-config_expected),real(driven(1)));

% Reject spin_config reference aliases on the finite-control dispatch as well
unsupported=supplied_config;
unsupported.reference=eye(2)/2;
rejected=false;
try
    controlled(model,2,unsupported,q,@finite_drive);
catch exception
    rejected=strcmp(exception.identifier,'cce_subsystem:configuration');
end
result=test_true(result,'Finite controls reject spin_config reference',rejected,...
                 'The driven evaluator must not ignore an unsupported normalization reference.');

% Compare integrated and Cartesian ideal rotations beside a finite drive
ideal=struct('spins',1,'time',0.15,'operator',0.3*z);
q=p;
q.refocusing={ideal};
accepted=true;
optional_error=Inf;
try
    implicit=controlled(model,2,spin_config,q,@finite_drive);
    q.refocusing={struct('spins',1,'time',0.15,'axis','z','angle',0.3)};
    empty_angle=controlled(model,2,spin_config,q,@finite_drive);
    q.refocusing={ideal};
    explicit=controlled(model,2,spin_config,q,@finite_drive);
    optional_error=max(norm(implicit-explicit),norm(empty_angle-explicit));
catch exception
    accepted=false;
    fprintf('operator representation exception %s\n',exception.identifier);
end
result=test_true(result,'Finite ideal operator optional angle',accepted,...
                 'Integrated and Cartesian native declarations preserve the same rotation.');
result=test_close(result,'Finite ideal operator angle parity',optional_error,0,1e-12,1e-12,...
                  'Equivalent native pulse declarations must match.');
fprintf('operator representation error %.3g\n',optional_error);

% Match integrated generators in Hamiltonian and master-equation paths
for method={'gcce','megcce'}
    q=rmfield(p,'drive');
    q.refocusing=ideal;
    if strcmp(method{1},'megcce')
        q.cce.basis.formalism='zeeman-liouv';
        q.bath_jumps=struct([]);
    end
    values=cell(1,3);
    for variant=1:3
        if variant==2
            q.refocusing=struct('spins',1,'time',0.15,'axis','z','angle',0.3);
        else
            q.refocusing=ideal;
        end
        values{variant}=clexpn(model.spin_system,@hahn_echo,q,'labframe');
    end
    result=test_close(result,[method{1} ' omitted ideal operator angle'],...
                      values{1},values{3},1e-12,1e-12,...
                      'An integrated Hermitian operator includes its physical angle.');
    result=test_close(result,[method{1} ' empty ideal operator angle'],...
                      values{2},values{3},1e-12,1e-12,...
                      'Cartesian and integrated native declarations agree.');
end

% Preserve normalised observable units through each generalised driver
for method={'gcce','megcce'}
    for driven=[false true]
        if driven&&strcmp(method{1},'megcce')
            continue;
        end
        q=p;
        q.output='observable';
        q.normalized=true;
        if strcmp(method{1},'megcce')
            q.cce.basis.formalism='zeeman-liouv';
            q.bath_jumps=struct([]);
        end
        q.coil={1e-20*eye(2),eye(2),1e20*eye(2)};
        sequence=@finite_drive;
        if ~driven
            q=rmfield(q,'drive');
            sequence=@free_evolution;
        end
        scaled=clexpn(model.spin_system,sequence,q,'labframe');
        result=test_close(result,[method{1} ' scaled observable driver ' num2str(driven)],...
                          scaled,ones(numel(q.times),3),1e-11,1e-11,...
                          'Operator units must not change normalized signal through any propagation path.');
    end
end

% Match disjoint one-site operator pulses across driven and ideal paths
for driven=[false true]
    q=p;
    ideal=struct('spins',2,'time',0.05,'operator',pi*x);
    if driven
        sequence=@finite_drive;
        q.refocusing={ideal};
    else
        sequence=@hahn_echo;
        q=rmfield(q,'drive');
        q.refocusing=ideal;
    end
    accepted=false;
    pulse_error=Inf;
    try
        [operator_result,operator_details]=clexpn(model.spin_system,sequence,q,'labframe');
        pulse=struct('spins',2,'time',0.05,'axis','x','angle',pi);
        if driven
            q.refocusing={pulse};
        else
            q.refocusing=pulse;
        end
        axis_result=clexpn(model.spin_system,sequence,q,'labframe');
        pulse_error=norm(operator_result(:)-axis_result(:));
        accepted=all(operator_details.reduction.valid,'all');
    catch exception
        fprintf('DISJOINT_PULSE driven=%d error=%s\n',driven,exception.identifier);
    end
    result=test_true(result,['Disjoint bath operator acceptance ' num2str(driven)],accepted,...
                     'A one-site exterior pulse has no effect when mean fields are disabled.');
    result=test_close(result,['Disjoint bath operator parity ' num2str(driven)],...
                      pulse_error,0,1e-12,1e-12,...
                      'Operator and axis-angle pulses must agree including the observed-spin-only baseline.');
    fprintf('DISJOINT_PULSE driven=%d accepted=%d error=%.3e\n',driven,accepted,pulse_error);
end

% Retain the explicit rejection when an exterior pulse can change mean fields
q=p;
q.refocusing={ideal};
q.cce.meanfield='static';
rejected=false;
try
    clexpn(model.spin_system,@finite_drive,q,'labframe');
catch exception
    rejected=strcmp(exception.identifier,'cce_segments:exterior_correlations');
end
result=test_true(result,'Exterior operator static-field rejection',rejected,...
                 'The no-mean-field skip must not discard a pulse that changes exterior expectations.');

% Require no-op control dispatch to preserve explicit and default normalisation
q=rmfield(p,'drive');
q.output='signal';
q.rho0=ones(2)/2;
ss=model.spin_system;
for normalised=[false true]
    q.normalized=normalised;
    plain=clexpn(ss,@free_evolution,q,'labframe');
    q.drive=struct('spins',1,'duration',0.01,'axis','x','amplitude',0,'detuning',0);
    driven=clexpn(ss,@finite_drive,q,'labframe');
    result=test_close(result,'Control dispatch preserves signal normalization',...
                      driven,plain,1e-12,1e-12,...
                      'A zero-amplitude control must not change the requested signal scale.');
    q=rmfield(q,'drive');
end
q.normalized=false;
plain=clexpn(ss,@free_evolution,q,'labframe');
q.drive=struct('spins',1,'duration',0.01,'axis','x','amplitude',0,'detuning',0);
driven=clexpn(ss,@finite_drive,q,'labframe');
result=test_close(result,'Control dispatch preserves default raw signal',...
                  driven,plain,1e-12,1e-12,...
                  'Omitting normalized must retain the same raw signal on both dispatch paths.');
fprintf('Control normalization dispatch discrepancy %.3g.\n',max(abs(driven-plain)));

% Compare native explicit and default Hahn declarations at the same signal scale
q=rmfield(q,'drive');
q.refocusing=struct('spins',1,'fraction',0.5,'operator',pi*x);
explicit_echo=clexpn(ss,@hahn_echo,q,'labframe');
q=rmfield(q,'refocusing');
q.cce.numerics.propagator='expm';
named_echo=clexpn(ss,@hahn_echo,q,'labframe');
result=test_close(result,'Named and explicit Hahn normalization',...
                  named_echo,explicit_echo,1e-12,1e-12,...
                  'Equivalent Hahn declarations must not change the requested signal scale.');

% Exercise the documented observable field for scalar and multiple detection operators
q.output='observable';
q.coil=full(z);
observed=clexpn(ss,@hahn_echo,q,'labframe');
q.output='density';
q=rmfield(q,'coil');
density_result=clexpn(ss,@hahn_echo,q,'labframe');
exact=zeros(numel(q.times),2);
for k=1:numel(q.times)
    exact(k,:)=[trace(z*reshape(density_result(k,:),2,2)) trace(x*reshape(density_result(k,:),2,2))];
end
result=test_close(result,'Driven scalar coil detection',observed,exact(:,1),...
                  1e-12,1e-12,'The documented operators field supplies the driven observable.');
q.output='observable';
q.coil={full(z),full(x)};
observed=clexpn(ss,@hahn_echo,q,'labframe');
result=test_close(result,'Driven cell coil detection',observed,exact,...
                  1e-12,1e-12,'Each observable has its own operator-norm error scale.');
fprintf('Driven operators discrepancy %.3g.\n',max(abs(observed-exact),[],'all'));

% Exercise the public driver through the named/finite dispatch hooks
integrated=clexpn(model.spin_system,@finite_drive,p,'labframe');
result=test_close(result,'Public finite-control driver',reshape(integrated(end,:),2,2),expected,...
                  1e-12,0,'Driver consumes the identical segment propagator');

% Approach the ideal pulse with a fixed rotation angle and decreasing width
ideal=struct('spins',1,'time',0.05,'axis','x','angle',pi);
p=rmfield(p,'drive');
p.refocusing={ideal};
[~,reference_aux]=controlled(model,2,spin_config,p,@refocusing_train);
errors=zeros(1,3);
p=rmfield(p,'refocusing');
for k=1:3
    width=10^(-k-2);
    p.drive=struct('spins',1,'time',0.05-width/2,...
                   'duration',width,'axis','x','amplitude',pi/width,'detuning',0);
    [~,s_aux]=controlled(model,2,spin_config,p,@finite_drive);
    errors(k)=norm(s_aux.density(:,:,end)-reference_aux.density(:,:,end),'fro');
end
result=test_true(result,'Finite pulse approaches ideal rotation',all(diff(errors)<0)&&errors(end)<2e-5,'Finite pulse approaches ideal rotation');

% Compare midpoint piecewise envelopes against a much finer independent product
mesh=[8 16 32];
waveform_errors=zeros(size(mesh));
mesh_count=1024;

% Accumulate midpoint propagators before tracing out the bath
spin=pauli(2);
x=spin.x;
U=eye(4);
dt=0.1/mesh_count;
for j=1:mesh_count
    U=expm(-1i*dt*(H+10*sin(pi*(j-0.5)/mesh_count)*kron(x,eye(2))))*U;
end
reference=partial(U*kron(p.rho0,spin_config.bath_rho{1})*U');
for k=1:numel(mesh)
    n=mesh(k);
    dt=0.1/n;
    p.drive=struct('spins',1,'time',(0:n-1)*dt,'duration',dt*ones(1,n),...
                   'axis','x','amplitude',10*sin(pi*((1:n)-0.5)/n),'detuning',0);
    [~,s_aux]=controlled(model,2,spin_config,p,@finite_drive);
    waveform_errors(k)=norm(s_aux.density(:,:,end)-reference,'fro');
end
result=test_true(result,'Piecewise midpoint refinement is second order',all(waveform_errors(1:end-1)./waveform_errors(2:end)>3.8),'Piecewise midpoint refinement is second order');

% Bound retained matrices before compiling a native XY8 train
q=rmfield(p,'drive');
q.times=linspace(0,0.1,100);
q.train=struct('name','xy8','count',8,'axis','x','spins',1);
q.cce.numerics.max_memory_bytes=1e5;
rejected=false;
try
    controlled(model,2,spin_config,q,@refocusing_train);
catch err
    rejected=strcmp(err.identifier,'cce_segments:memory');
end
result=test_true(result,'Retained control segment memory cap',rejected,...
                 'The memory cap must cover all retained Hamiltonians and pulses across requested times.');

% Apply the same retained-control bound through the public scheduled driver
q.cce.numerics.max_memory_bytes=5e5;
rejected=false;
try
    clexpn(model.spin_system,@refocusing_train,q,'labframe');
catch err
    rejected=strcmp(err.identifier,'cce_schedule:memory');
end
result=test_true(result,'Scheduled retained control memory cap',rejected,...
                 'The public driver must reserve the time-by-segment storage on active workers.');

% Check literal XY8 phases and analytic Uhrig timestamps
p=rmfield(p,'drive');
p.train=struct('name','xy8','count',8,'axis','x','spins',1);
[~,s_aux]=controlled(model,2,spin_config,p,@refocusing_train);
result=test_close(result,'XY8 timestamps',s_aux.segments.schedule.positions(:,end),0.1*((1:8)'-0.5)/8,1e-15,0,'XY8 timestamps');
result=test_true(result,'XY8 phases',strcmp([s_aux.segments.schedule.events.axis],'xyxyyxyx'),'XY8 phases');
p.train=struct('name','udd','count',5,'axis','y','spins',1);
[~,s_aux]=controlled(model,2,spin_config,p,@refocusing_train);
result=test_close(result,'UDD timestamps',s_aux.segments.schedule.positions(:,end),0.1*sin(pi*(1:5)'/12).^2,1e-15,0,'UDD timestamps');

% Exercise time-dependent exterior mean fields and explicit finite-drive refusal
p=rmfield(p,'train');
p.cce.meanfield='static';
p.refocusing={ideal};
p.refocusing{1}.spins=2;
[~,s_aux]=controlled(model,[],spin_config,p,@refocusing_train);
field=s_aux.segments.exterior_expectations{end};
result=test_close(result,'Exterior ideal pulse flips the sampled field',field{2}(1,3),-spin_config.expectations(1,3),1e-14,0,'Exterior ideal pulse flips the sampled field');
h_before=s_aux.segments.hamiltonians{end}{2};
h_after=s_aux.segments.hamiltonians{end}{3};
result=test_true(result,'Exterior pulse changes subsequent Hamiltonian',norm(h_before-h_after,'fro')>1,'Exterior pulse changes subsequent Hamiltonian');
p=rmfield(p,'refocusing');
p.drive=struct('spins',2,'duration',0.1,'axis','x','amplitude',7,'detuning',0);
refused=false;
try
    controlled(model,[],spin_config,p,@finite_drive);
catch exception
    refused=strcmp(exception.identifier,'finite_drive:spins');
end
result=test_true(result,'Finite exterior drive is refused rather than frozen',refused,'Finite exterior drive is refused rather than frozen');

% Refuse finite observed controls in the projected-manifold approximation
p.cce.meanfield='none';
p.cce.approximation='projected';
p.drive.spins=1;
refused=false;
try
    clexpn(model.spin_system,@finite_drive,p,'labframe');
catch exception
    refused=strcmp(exception.identifier,'clexpn:projection');
end
result=test_true(result,'Projected finite observed drive refused',refused,'General observed drives need gCCE');

% Verify an explicitly ordered joint pulse against a dense two-spin operator
p.cce.approximation='generalized';
p=rmfield(p,'drive');
p.refocusing={struct('spins',[2 1],'time',0,'operator',0.7*kron(z,x))};
[~,s_aux]=controlled(model,2,spin_config,p,@refocusing_train);
expected_pulse=expm(-1i*0.7*kron(x,z));
result=test_close(result,'Selective observed-spin system operator ordering',s_aux.segments.pulses{end}{1},expected_pulse,...
                  1e-13,0,'The target tensor order is remapped to canonical physical order');

% Distinguish a physical spin-one rotation from a selective two-level pulse
sys.isotopes={'14N'};
sys.magnet=0;
sys.disable={'hygiene'};
sys.parallel={'processes',1};
sys.output='hush';
inter=struct();
spin_system=create(sys,inter);
q=rmfield(p,'transition');
q.observed=1;
q.cce.bath_spins=[];
q.rho0=diag([1 0 0]);
q.cce.order=0;
q.cce.bath_state=struct('mode','mixed_product');
q.refocusing={struct('spins',1,'time',0,'axis','x','angle',pi)};
q=cce_options(spin_system,q,@refocusing_train,'labframe');
spin_one=cce_compile(spin_system,q,'labframe');
states=struct('bath_ids',[],'bath_rho',{{}});
[~,physical_aux]=controlled(spin_one,[],states,q,@refocusing_train);
q.refocusing={struct('spins',1,'time',0,'operator',pi*[0 1 0;1 0 0;0 0 0]/2)};
[~,selective_aux]=controlled(spin_one,[],states,q,@refocusing_train);
result=test_close(result,'Physical spin-one pi pulse',physical_aux.density(:,:,end),diag([0 0 1]),...
                  1e-12,0,'A full spin-one rotation maps m=+1 to m=-1');
result=test_close(result,'Selective two-level pi pulse',selective_aux.density(:,:,end),diag([0 1 0]),...
                  1e-12,0,'The declared two-level pulse maps m=+1 to m=0');
result=test_true(result,'Physical and selective pulses differ',...
                 norm(physical_aux.density(:,:,end)-selective_aux.density(:,:,end),'fro')>1,...
                 'A spin-one Cartesian rotation must not silently become a two-level pulse');

% Apply the declared integrated generator on the static native path
sys=struct();
sys.magnet=0;
sys.isotopes={'E'};
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter=struct();
inter.zeeman.matrix={zeros(3)};
spin_system=create(sys,inter);
stat.observed=1;
stat.cce.bath_spins=[];
stat.cce.approximation='generalized';
stat.cce.expansion='standard';
stat.output='density';
stat.cce.order=0;
stat.rho0=diag([1 0]);
stat.transition=struct('spins',1,'ms',[0.5 -0.5]);
stat.orientation=[0 0 0];
stat.cce.interaction_cutoff=Inf;
stat.cce.graph=struct('mode','complete');
stat.cce.meanfield='none';
stat.times=[0 1];
stat.cce.bath_state=struct('mode','mixed_product');
stat.cce.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
stat.refocusing=struct('spins',1,'fraction',0.5,'operator',pi*[0 1;1 0]/2);
rotated=clexpn(spin_system,@hahn_echo,stat,'labframe');
stat.refocusing.operator=zeros(2);
unrotated=clexpn(spin_system,@hahn_echo,stat,'labframe');
stat.refocusing=struct('spins',1,'fraction',0.5,'axis','x','angle',pi);
folded=clexpn(spin_system,@hahn_echo,stat,'labframe');
result=test_close(result,'Static explicit generator with angle pi',reshape(rotated(end,:),2,2),diag([0 1]),...
                  1e-12,0,'An integrated pi*Sx generator inverts the observed population');
result=test_close(result,'Static explicit generator with zero angle',reshape(unrotated(end,:),2,2),diag([1 0]),...
                  1e-12,0,'A zero generator must leave the observed state unchanged');
result=test_close(result,'Static angle-folded generator',reshape(folded(end,:),2,2),reshape(rotated(end,:),2,2),...
                  1e-12,0,'Integrated pi*Sx equals a Cartesian pulse with angle pi');

% Print numerical evidence without adding nonstandard harness fields
fprintf('controls dense error %.3g; ideal-limit errors %s; waveform errors %s\n',error,mat2str(errors,6),mat2str(waveform_errors,6));

% Pulse-generator validation is unchanged by reciprocal operator-angle scaling
[model,q,ensemble]=fixture();
spin_config=ensemble.spin_configs(1);
rejected=false(1,2);
q=rmfield(q,'drive');
q.refocusing={struct('spins',1,'fraction',0.5,'operator',1e20*(1e-20*[0 1;0 0]))};
try
    controlled(model,2,spin_config,q,@refocusing_train);
catch
    rejected(1)=true;
end
try
    clexpn(model.spin_system,@refocusing_train,q,'labframe');
catch
    rejected(2)=true;
end
result=test_true(result,'Scaled non-Hermitian pulse rejected',all(rejected),...
                 'A small operator norm cannot hide non-Hermiticity compensated by a large angle.');
fprintf('R6 scaled pulse rejections %d/2\n',sum(rejected));

% Preserve transition-derived detection through both segmented control routes
for coil_sequence={@refocusing_train,@finite_drive}
    [coil_model,coil_request,~]=fixture();
    coil_request.output='observable';
    coil_request=rmfield(coil_request,'drive');
    coil_request.refocusing={struct('spins',1,'time',0.05,'axis','x','angle',pi)};
    if isequal(coil_sequence{1},@finite_drive)
        coil_request.drive=struct('spins',1,'duration',0.1,'axis','x','amplitude',7,'detuning',3);
    end
    implicit=clexpn(coil_model.spin_system,coil_sequence{1},coil_request,'labframe');
    coil_request.coil=coil_model.transition.alpha*coil_model.transition.beta';
    explicit=clexpn(coil_model.spin_system,coil_sequence{1},coil_request,'labframe');
    result=test_close(result,['Optional coil: ' func2str(coil_sequence{1})],implicit,explicit,...
                      1e-12,1e-12,'Omitted coils use the selected transition in segmented propagation.');
end

end

% Trace the independent bath index in the dense Kronecker reference
function rho=partial(full)

% Sum the bath diagonal blocks
rho=zeros(2);
for a=1:2
    for b=1:2
        rho(a,b)=trace(full((a-1)*2+(1:2),(b-1)*2+(1:2)));
    end
end

end

% Build a small noncommuting observed-spin-bath model in angular-frequency units
function [model,p,ensemble]=fixture()

% Compile the explicit two-spin model and its product state
sys.isotopes={'1H','1H'};
sys.magnet=0;
sys.disable={'hygiene'};
sys.parallel={'processes',1};
sys.output='hush';
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([2 3 5])/(2*pi);
spin_system=create(sys,inter);
p.observed=1;
p.cce.bath_spins=2;
p.cce.approximation='generalized';
p.cce.expansion='standard';
p.output='density';
p.times=[0 0.1];
p.cce.order=1;
p.rho0=[1 1;1 1]/2;
p.orientation=[0 0 0];
p.normalized=false;
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','complete','score','coupling_norm');
p.cce.bath_state=struct('mode','explicit_product','rho',{{diag([1 0])}});
p.cce.meanfield='none';
p.cce.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
p.drive=struct('spins',1,'duration',0.1,'axis','x','amplitude',7,'detuning',3);
p=cce_options(spin_system,p,@finite_drive,'labframe');
model=cce_compile(spin_system,p,'labframe');
configuration=struct('bath_ids',2,'bath_rho',{{diag([1 0])}},'state_id',1,'configuration_id',1,'weight',1);
local=cce_subsystem(model,2,configuration,[]);
ensemble.spin_configs=rmfield(local.configuration,{'exterior_ids','exterior_expectations'});

end

% Evaluate native experiments directly on an explicitly conditioned subsystem
function [answer,aux]=controlled(model,sites,configuration,p,experiment)

% Construct the conditioned subsystem and execute its native experiment
p=cce_options(model.spin_system,p,experiment,'labframe');
model=cce_compile(model.spin_system,p,'labframe');
local=cce_subsystem(model,sites,configuration,[]);
[H,R,K,p]=cce_generators(local,p,[]);
[answer,aux]=experiment(local,p,H,R,K);

end


