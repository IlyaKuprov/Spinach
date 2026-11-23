% Test: kernel/clexpn_master | Native master equation CCE
% Native master-equation, projected-block, and conditioning callback tests.
% Syntax:
%
%                    result=test_clexpn_master()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_master()

% Specify the physical system and regression targets
result=new_test_result('kernel/clexpn_master','Native master equation CCE',...
                       'GKSL evolution and conditioned evaluator preservation.');
sys.magnet=0;
sys.isotopes={'E','1H','1H'};
sys.output='hush';
sys.parallel={'processes',1};
sys.tols.cce_atol=1e-10;
sys.tols.cce_rtol=1e-10;
sys.tols.cce_mem=2^26;
pool=gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
inter.zeeman.scalar={0,0,0};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 0.7]);
inter.coupling.matrix{1,3}=diag([0 0 0.3]);
s=create(sys,inter);
p.observed=1;
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.times=[0 0.1 0.2];
p.bath_jumps=struct('support',2,'operator','L-','rate',0.4);
p.cce.order=2;
[a,da]=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Complete dissipative bath',a,da.raw{1}(end).value,1e-12,1e-12,...
                  'Full-order reduction agrees with direct raw GKSL propagation.');
p.cce.approximation='projected';
b=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Projected coherence block',a,b,1e-12,1e-12,...
                  'Bath-only relaxation preserves the selected observed manifolds.');
projected=clexpn(s,@hahn_echo,p,'labframe');
p.cce.approximation='generalized';
generalized=clexpn(s,@hahn_echo,p,'labframe');
result=test_close(result,'Dissipative Hahn echo',projected,generalized,1e-12,1e-12,...
                  'Both representations use the same fixed physical receiver.');
echo_parameters=p;
echo_parameters.output='density';
[echo_density,echo_details]=clexpn(s,@hahn_echo,echo_parameters,'labframe');
result=test_close(result,'Hahn density full bath',echo_density,echo_details.raw{1}(end).value,...
                  1e-12,1e-12,'Generalized echo supports direct observed-density output.');
p.cce.expansion='standard';
p.cce.order=1;
p.bath_samples=2;
p.bath_seed=7;
p.cce.interlaced=struct('max_states',4);
[conditional,di]=clexpn(s,@free_evolution,p,'labframe');
result=test_true(result,'Conditioned factors finite',all(isfinite(conditional)),...
                 'Recursive factors use one common parent conditioning event.');
result=test_true(result,'Irreducible metadata',strcmp(di.raw{1}(2).kind,'irreducible'),...
                 'The reducer consumes an irreducible factor exactly once.');
p.cce.order=2;
[complete,df]=clexpn(s,@free_evolution,p,'labframe');
p.cce.expansion='standard';
p=rmfield(p,{'bath_samples','bath_seed'});
p.cce=rmfield(p.cce,'interlaced');
exact=clexpn(s,@free_evolution,p,'labframe');
result=test_close(result,'Complete conditioned bath',complete,exact,1e-12,1e-12,...
                  'Full bath conditioning reaches the supplied dissipative evaluator.');
result=test_true(result,'Full enumeration',df.raw{1}(end).assignment_count==4,...
                 'All physical assignments enter the complete-bath callback.');
result=test_true(result,'Full-bath lower orders not evaluated',...
                 ~any(df.reduction.by_order_valid(:,:,1:end-1),'all'),...
                 'The exact complete-bath shortcut does not report artificial lower-order convergence.');

% Observed-spin relaxation is part of the retained physical GKSL contract
relax=p;
relax.cce.bath_spins=[];
relax.cce.order=0;
relax.cce.expansion='standard';
relax.cce.approximation='generalized';
relax.output='density';
relax.bath_jumps=struct('support',1,'operator','L+','rate',0.4);
rho=clexpn(s,@free_evolution,relax,'labframe');
result=test_close(result,'Observed relaxation population',rho(:,4),0.5*exp(-0.4*p.times(:)),...
                  1e-12,1e-12,'Observed-spin relaxation is not silently discarded.');

% Zero-rate manifold-changing jumps vanish before projected closure checks
zero=relax;
zero.cce.approximation='projected';
zero.output='signal';
zero.bath_jumps=struct('support',1,'operator',[0 1;0 0],'rate',0);
zero_model_p=cce_options(s,zero,@cce_compile,'labframe');
zero_model=cce_compile(s,zero_model_p,'labframe');
[zero_system,~]=cce_subsystem(zero_model,zero_model_p,[],struct(),[]);
zero_parameters=zero_model_p;
[with_zero,zero_loss,~]=cce_generators(zero_system,zero_parameters,[]);
zero_parameters.bath_jumps=struct('support',{},'operator',{},'rate',{});
[without_jump,~,~]=cce_generators(zero_system,zero_parameters,[]);
zero_generator=-1i*(with_zero{1}-with_zero{2})+zero_loss;
no_generator=-1i*(without_jump{1}-without_jump{2});
result=test_close(result,'Zero-rate projected jump omitted',zero_generator,no_generator,0,0,...
                  'A zero coefficient contributes neither dissipation nor manifold leakage.');

% Active exterior jumps cannot supply a frozen mean field
exterior=zero;
exterior.cce.bath_spins=2;
exterior.cce.order=1;
exterior.cce.meanfield='static';
exterior.bath_jumps.support=2;
exterior.bath_jumps.rate=1;
exterior_model_p=cce_options(s,exterior,@cce_compile,'labframe');
exterior_model=cce_compile(s,exterior_model_p,'labframe');
configuration=struct('bath_ids',2,'bath_rho',{{diag([1 0])}},...
                     'exterior_ids',2,'exterior_expectations',[0 0 0.5]);
[exterior_system,~]=cce_subsystem(exterior_model,exterior_model_p,[],configuration,[]);
rejected=false;
try
    cce_generators(exterior_system,exterior_model_p,[]);
catch exception
    rejected=strcmp(exception.identifier,'cce:exteriorJump');
end
result=test_true(result,'Active exterior jump refused',rejected,...
                 'An evolving exterior population cannot silently become a static field.');

% Exercise the recursive and complete-bath conditioned paths
model_p=cce_options(s,df.parameters,@cce_compile,'labframe');
model=cce_compile(s,model_p,'labframe');
model.bas.formalism='zeeman-liouv';
model_clusters=df.clusters;
local=df.parameters;
local.cce.interlaced.reference_state=1;
raw=cce_interlace(model,df.clusters(1),model_clusters,df.ensemble,@free_evolution,local);
result=test_true(result,'Recursive conditioned factor',strcmp(raw.kind,'irreducible'),...
                 'The empty reference and denominator subsets enter the irreducible factor.');
raw=cce_interlace(model,df.clusters(end),model_clusters,df.ensemble,@free_evolution,local);
result=test_true(result,'Full-bath conditioned factor',raw.full_bath&&raw.assignment_count==4,...
                 'The complete-bath shortcut evaluates every physical assignment.');

% Checkpoint payload retains contract version, moments, masks, and uncertainties
file=[tempname '.mat'];
p.cce.checkpoint=struct('file',file,'resume',false);
[first,~]=clexpn(s,@free_evolution,p,'labframe');
p.cce.checkpoint.resume=true;
[resumed,~]=clexpn(s,@free_evolution,p,'labframe');
delete(file);
result=test_close(result,'Native checkpoint resume',first,resumed,0,0,'Exact task records survive restart.');
p.cce=rmfield(p.cce,'checkpoint');
delete(gcp('nocreate'));
sys.parallel={'processes',2};
parallel_system=create(sys,inter);
parallel=clexpn(parallel_system,@free_evolution,p,'labframe');
result=test_close(result,'Native parallel determinism',parallel,first,0,0,'Workers return raw indexed records only.');

% Finite Liouville drives retain uncertainty under strict numerical budgets
request=p;
delete(gcp('nocreate'));
sys.parallel={'processes',1};
strict=create(sys,inter);
strict.tols.cce_atol=1e-30;
strict.tols.cce_rtol=0;
request.times=[0 0.05 0.1 0.2];
request.drive=struct('spins',1,'axis','x','time',0.02,...
                     'duration',0.1,'amplitude',7,'detuning',2);
request.drive(2)=request.drive;
request.drive(2).time=0.06;
request.drive(2).axis='y';
request.drive(2).amplitude=3;
[~,driven]=clexpn(strict,@finite_drive,request,'labframe');
request=driven.parameters;
drive_model_p=cce_options(strict,request,@cce_compile,'labframe');
drive_model=cce_compile(strict,drive_model_p,'labframe');
drive_model.bas.formalism='zeeman-liouv';
checks=cce_diagnostics(drive_model,driven.raw,driven.reduction,request);
result=test_true(result,'Driven Liouville numerical budget is not exact',...
                 all(driven.raw{1}(end).abs_error_bound(2:end)>0)&&...
                 ~any(checks.numerical_budget_met(2:end),'all'),...
                 'Full-bath reductions propagate nonzero segment estimates into budget diagnostics.');

% Independently sum the declared arithmetic estimate over all finite intervals
[retained,~]=cce_subsystem(drive_model,drive_model_p,drive_model_p.cce.bath_spins,struct(),[]);
[H,R,K,retained]=cce_generators(retained,request,[]);
local=request;
L=H+1i*R+1i*K;
[~,aux]=finite_drive(retained,local,H,R,K);
hilbert=retained;
hilbert.bas.formalism='zeeman-hilb';
x=cce_operator(hilbert,'Lx',1);
y=cce_operator(hilbert,'Ly',1);
z=cce_operator(hilbert,'Lz',1);
identity=speye(retained.cce.dimension);
controls={7*x+2*z,3*y+2*z};
for k=1:2
    controls{k}=kron(identity,controls{k})-kron(controls{k}.',identity);
end
estimates=zeros(numel(request.times),1);
for n=1:numel(request.times)
    time=request.times(n);
    starts=[request.drive.time];
    stops=starts+[request.drive.duration];
    cuts=unique([0 time starts(starts>0&starts<time) stops(stops>0&stops<time)]);
    for k=1:numel(cuts)-1
        generator=L;
        midpoint=(cuts(k)+cuts(k+1))/2;
        for j=find(starts<=midpoint&stops>midpoint)
            generator=generator+controls{j};
        end
        estimates(n)=estimates(n)+1+(cuts(k+1)-cuts(k))*norm(generator,1);
    end
end
estimates=64*eps*size(L,1)*estimates;
result=test_close(result,'Driven Liouville segment estimate',aux.abs_error_bound,2*estimates,...
                  16*eps*norm(estimates),0,...
                  'Selected normalised coherence divides by its one-half preparation reference.');

% Explicit receiver magnitude and normalisation scale the same propagation estimate
local.coil={[0 3;0 0];[0 0;0 2]};
local.normalized=false;
[~,receivers]=finite_drive(retained,local,H,R,K);
result=test_close(result,'Driven Liouville receiver scaling',receivers.abs_error_bound,...
                  estimates*[3 2],8*eps*norm(estimates*[3 2],'fro'),0,...
                  'Observable bounds follow receiver-column order.');
local.normalized=true;
[~,normalised]=finite_drive(retained,local,H,R,K);
result=test_close(result,'Driven Liouville reference scaling',normalised.abs_error_bound,...
                  estimates*[2 2],8*eps*norm(estimates*[2 2],'fro'),0,...
                  'Each receiver divides by its own initial expectation.');
local.output='density';
[~,densities]=finite_drive(retained,local,H,R,K);
result=test_close(result,'Driven Liouville density scaling',densities.abs_error_bound,...
                  repmat(estimates,1,4),16*eps*norm(estimates),0,...
                  'Density entries retain the unscaled segment estimate.');

end


