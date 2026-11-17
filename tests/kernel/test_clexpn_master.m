% Native master-equation, projected-block and conditioning callback tests.
% Syntax: result=test_clexpn_master()
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_master()
result=new_test_result('kernel/clexpn_master','Native master equation CCE',...
                       'GKSL evolution and conditioned evaluator preservation.');
sys.magnet=0; sys.isotopes={'E','1H','1H'}; sys.output='hush'; sys.parallel={'local',2};
inter.zeeman.scalar={0,0,0}; inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 0.7]); inter.coupling.matrix{1,3}=diag([0 0 0.3]);
s=create(sys,inter); p.observed=1;
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]); p.times=[0 0.1 0.2];
p.bath_jumps=struct('support',2,'operator','L-','rate',0.4);
p.cce.order=2; p.cce.basis=struct('formalism','zeeman-liouv');
p.cce.numerics=struct('atol',1e-10,'rtol',1e-10,'max_memory_bytes',2^26);
[a,da]=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Complete dissipative bath',a,da.raw{1}(end).value,1e-12,1e-12,...
                  'Full-order reduction agrees with direct raw GKSL propagation.');
p.cce.approximation='projected'; b=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Projected coherence block',a,b,1e-12,1e-12,...
                  'Bath-only relaxation preserves the selected observed manifolds.');
projected=clexpn(s,@hahn_echo,p,'labframe'); p.cce.approximation='generalized';
generalized=clexpn(s,@hahn_echo,p,'labframe');
result=test_close(result,'Dissipative Hahn echo',conj(projected),generalized,1e-12,1e-12,...
                  'The swapped projected block is the conjugate coherence selected by the physical receiver.');
echo_parameters=p; echo_parameters.output='density';
[echo_density,echo_details]=clexpn(s,@hahn_echo,echo_parameters,'labframe');
result=test_close(result,'Hahn density full bath',echo_density,echo_details.raw{1}(end).value,...
                  1e-12,1e-12,'Generalized echo supports direct observed-density output.');
p.cce.expansion='standard'; p.cce.order=1;
p.cce.bath_state=struct('mode','interlaced','basis','uniform','probabilities','uniform',...
                        'samples',2,'seed',7); p.cce.interlaced=struct('max_states',4);
[conditional,di]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Conditioned factors finite',all(isfinite(conditional)),...
                 'Recursive factors use one common parent conditioning event.');
result=test_true(result,'Irreducible metadata',strcmp(di.raw{1}(2).kind,'irreducible'),...
                 'The reducer consumes an irreducible factor exactly once.');
p.cce.order=2; [complete,df]=clexpn(s,@free_evolution,p,'labframe');
p.cce.expansion='standard'; p.cce.bath_state=struct('mode','mixed_product');
p.cce=rmfield(p.cce,'interlaced'); exact=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Complete conditioned bath',complete,exact,1e-12,1e-12,...
                  'Full bath conditioning reaches the supplied dissipative evaluator.');
result=test_true(result,'Full enumeration',df.raw{1}(end).assignment_count==4,...
                 'All physical assignments enter the complete-bath callback.');
result=test_true(result,'Full-bath lower orders not evaluated',...
                 ~any(df.reduction.by_order_valid(:,:,1:end-1),'all'),...
                 'The exact complete-bath shortcut does not report artificial lower-order convergence.');

% Observed-spin relaxation is part of the retained physical GKSL contract
relax=p; relax.cce.bath_spins=[]; relax.cce.order=0;
relax.cce.expansion='standard'; relax.cce.bath_state=struct('mode','mixed_product');
if isfield(relax.cce,'interlaced'), relax.cce=rmfield(relax.cce,'interlaced'); end
relax.cce.approximation='generalized';
relax.output='density'; relax.bath_jumps=struct('support',1,'operator','L+','rate',0.4);
rho=clexpn(s,@free_evolution,relax,'labframe');
result=test_close(result,'Observed relaxation population',rho(:,4),0.5*exp(-0.4*p.times(:)),...
                  1e-12,1e-12,'Observed-spin relaxation is not silently discarded.');

% Zero-rate manifold-changing jumps vanish before projected closure checks
zero=relax; zero.cce.approximation='projected'; zero.output='signal';
zero.bath_jumps=struct('support',1,'operator',[0 1;0 0],'rate',0);
zero_model=cce_compile(s,zero,'labframe');
[zero_system,~]=cce_subsystem(zero_model,[],struct(),[]);
zero_parameters=zero_model.parameters;
[with_zero,~,~]=cce_generators(zero_system,zero_parameters,[]);
zero_parameters.bath_jumps=struct('support',{},'operator',{},'rate',{});
[without_jump,~,~]=cce_generators(zero_system,zero_parameters,[]);
zero_generator=-1i*(with_zero.left-with_zero.right)+with_zero.dissipator;
no_generator=-1i*(without_jump.left-without_jump.right);
result=test_close(result,'Zero-rate projected jump omitted',zero_generator,no_generator,0,0,...
                  'A zero coefficient contributes neither dissipation nor manifold leakage.');

% Active exterior jumps cannot supply a frozen mean field
exterior=zero; exterior.cce.bath_spins=2; exterior.cce.order=1;
exterior.cce.meanfield='static'; exterior.bath_jumps.support=2;
exterior.bath_jumps.rate=1; exterior_model=cce_compile(s,exterior,'labframe');
configuration=struct('bath_ids',2,'bath_rho',{{diag([1 0])}},...
                     'exterior_ids',2,'exterior_expectations',[0 0 0.5]);
[exterior_system,~]=cce_subsystem(exterior_model,[],configuration,[]);
rejected=false;
try
    cce_generators(exterior_system,exterior_model.parameters,[]);
catch exception
    rejected=strcmp(exception.identifier,'cce:exteriorJump');
end
result=test_true(result,'Active exterior jump refused',rejected,...
                 'An evolving exterior population cannot silently become a static field.');

% Instrument the supplied evaluator on recursive and complete-bath paths
model=cce_compile(s,df.parameters,'labframe'); model.formalism='zeeman-liouv';
model.clusters=df.clusters; calls=[]; events={}; local=df.parameters;
local.cce.interlaced.reference_state=1; schedule.pulse_sequence=@free_evolution;
schedule.contract=cce_capabilities(@free_evolution);
callback=schedule.contract.evaluator;
cce_interlace(model,df.clusters(1),df.ensemble,schedule,local,@counted);
result=test_true(result,'Recursive evaluator callback',any(calls==0)&&any(calls==df.clusters(1).id),...
                 'The empty reference and denominator subsets invoke the supplied evaluator.');
calls=[]; cce_interlace(model,df.clusters(end),df.ensemble,schedule,local,@counted);
result=test_true(result,'Full-bath evaluator callback',numel(calls)==4&&all(calls==df.clusters(end).id),...
                 'The complete-bath shortcut never substitutes a closed-system evaluator.');
result=test_true(result,'Conditioned physical assignments',~isempty(events),...
                 'Callbacks receive explicit parent-conditioned bath densities.');


% Checkpoint payload retains contract version, moments, masks and uncertainties
file=[tempname '.mat']; p.cce.checkpoint=struct('file',file,'resume',false);
[first,~]=clexpn(s,@free_evolution,p,'labframe'); p.cce.checkpoint.resume=true;
[resumed,~]=clexpn(s,@free_evolution,p,'labframe'); delete(file);
result=test_close(result,'Native checkpoint resume',first,resumed,0,0,'Exact task records survive restart.');
p.cce=rmfield(p.cce,'checkpoint'); p.cce.parallel=struct('mode','clusters','batch_size',2);
parallel=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Native parallel determinism',parallel,first,0,0,'Workers return raw indexed records only.');

% Capture physical conditioning, not just callback counts
function raw=counted(model,cluster,configuration,schedule,parameters)
calls(end+1)=cluster.id; events{end+1}=configuration.bath_rho;
raw=callback(model,cluster,configuration,schedule,parameters);
end
end


