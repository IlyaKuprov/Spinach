% Density physicality and singular-baseline exact bypass tests. Syntax:
%
%                         result=test_gcce_density()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_density()

% State the no-repair requirement for elementwise density expansions
result=new_test_result('kernel/gcce_density','gCCE density diagnostics',...
                       'Baseline zeros mask truncated ratios and full order bypasses them.');
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
spin_system=create(sys,inter);
spin_system.inter.zeeman.matrix{1}(3,1)=100;
spin_system.inter.coupling.matrix{1,2}=diag([31 43 57]);
p.observed=1;
p.cce.bath_spins=2;
p.cce.approximation='generalized';
p.cce.expansion='standard';
p.output='density';
p.rho0=[1 0;0 0];
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','complete','score','coupling_norm');
p.cce.meanfield='none';
p.times=[0 pi/200 pi/100];
p.cce.bath_state=struct('mode','mixed_product');
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Exact full-bath density must stay finite at zero observed-spin-only elements
p.cce.order=1;
[actual,actual_diagnostics]=clexpn(spin_system,@free_evolution,p,'labframe');
density=reshape(actual.',2,2,[]);
trace_error=zeros(numel(p.times),1);
hermiticity_error=trace_error;
minimum_eigenvalue=trace_error;
for k=1:numel(p.times)
    rho=density(:,:,k);
    trace_error(k)=abs(trace(rho)-1);
    hermiticity_error(k)=norm(rho-rho','fro');
    minimum_eigenvalue(k)=min(eig((rho+rho')/2));
end
result=test_true(result,'full-order bypass finite',all(actual_diagnostics.reduction.valid,'all'),...
                 'Exact reduced density does not require a baseline division.');
result=test_close(result,'initial density',density(:,:,1),p.rho0,1e-12,0,...
                  'Zero duration returns the supplied observed density.');
result=test_true(result,'unit trace',max(trace_error)<1e-12,...
                 'Unitary full-system propagation preserves partial-trace normalization.');
result=test_true(result,'Hermitian',max(hermiticity_error)<1e-12,...
                 'No Hermitian repair is applied.');
result=test_true(result,'positive',min(minimum_eigenvalue)>-1e-12,...
                 'Exact partial trace is positive without eigenvalue clipping.');

% Force a truncated family with a zero baseline and retain the raw record
p.cce.bath_spins=2;
p.cce.order=0;
truncated=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_true(result,'observed baseline preserved',all(isfinite(truncated),'all'),...
                 'Order zero returns the direct observed-spin-only density.');
p=actual_diagnostics.parameters;
clusters=actual_diagnostics.clusters;
records=actual_diagnostics.raw;
ratio_parameters=p;
ratio_parameters.cce.order=2;
ratio_parameters.cce.bath_spins=[2 3];
ratio_reduced=cce_reduce(records,clusters,ratio_parameters);
result=test_true(result,'elementwise baseline zeros masked',any(~ratio_reduced.valid,'all'),...
                 'The reducer alone must not replace undefined ratios by one.');
result=test_true(result,'raw density retained',all(isfinite(records{1}(2).value),'all'),...
                 'Masking factors does not destroy direct raw density information.');

% Select one bath spin but retain the polarised exterior in the static field
sys.isotopes{3}='1H';
inter.zeeman.matrix{3}=zeros(3);
inter.coupling.matrix=cell(3);
spatial=create(sys,inter);
spatial.inter.coordinates={[0 0 0],[1 0 0],[10 0 0]};
spatial.inter.coupling.matrix{1,2}=diag([4 0 2]);
spatial.inter.coupling.matrix{1,3}=[0 0 0;0 0 0;14 0 0];
radius=p;
radius.cce.bath_spins=[2 3];
radius.cce.bath_radius=2;
radius.cce.bath_state=struct('mode','explicit_product','rho',{{eye(2)/2,diag([1 0])}});
radius.cce.meanfield='static';
[up,retained]=clexpn(spatial,@free_evolution,radius,'labframe');
result=test_close(result,'Complete radius-selected density',up,retained.raw{1}(2).value,0,0,...
                  'Every element, including exact initial zeros, comes from the retained cluster.');
radius.cce.bath_state.rho{2}=diag([0 1]);
down=clexpn(spatial,@free_evolution,radius,'labframe');
result=test_true(result,'Radius-excluded exterior field',imag(up(end,2))*imag(down(end,2))<0,...
                 'Reversing only the excluded spin reverses the observed transverse coherence.');
radius.cce.bath_radius=10;
[~,incomplete]=clexpn(spatial,@free_evolution,radius,'labframe');
result=test_true(result,'Incomplete radius-selected density',...
                 isequal(incomplete.graph.bath_ids,[2 3])&&~any(incomplete.reduction.valid(1,2:4)),...
                 'Two retained sites with no joint cluster still have undefined zero-baseline ratios.');

% Reject nonpositive densities and unresolved initial scalar normalisation
bad=p;
bad.rho0=diag([1.1 -0.1]);
caught=false;
try
    cce_compile(spin_system,bad,'labframe');
catch
    caught=true;
end
result=test_true(result,'negative input rejected',caught,'State positivity is checked.');
bad=p;
bad.output='signal';
bad.normalized=true;
bad.transition.spins=1;
bad.transition.alpha=[1;0];
bad.transition.beta=[0;1];
caught=false;
try
    clexpn(spin_system,@free_evolution,bad,'labframe');
catch
    caught=true;
end
result=test_true(result,'zero normalization rejected',caught,...
                 'Initial zero signal has no normalized scalar interpretation.');

% Reproduce fixed and sampled product states without changing the global stream
sample=p;
sample.cce.bath_state=struct('mode','sampled','samples',5,'seed',17,...
                             'basis','uniform','probabilities',{{[1 0]}});
sample.cce.averaging='after_expansion';
before=rng;
[~,ensemble]=clexpn(spin_system,@free_evolution,sample,'labframe');
after=rng;
[~,again]=clexpn(spin_system,@free_evolution,sample,'labframe');
result=test_true(result,'global RNG unchanged',isequal(before,after),...
                 'Sampling owns counter-based Threefry substreams.');
result=test_true(result,'fixed probability state',all(cellfun(@(rho)isequal(rho,diag([1 0])),...
                     [ensemble.ensemble.spin_config.bath_rho])),...
                 'Probability-one states are not resampled into forbidden levels.');
result=test_true(result,'reproducible assignments',isequal(ensemble.ensemble.spin_config,again.ensemble.spin_config),...
                 'Explicit state IDs determine the stream independently of workers.');
sample.cce.bath_state.probabilities={[0.5 0.5]};
sample.cce.meanfield='static';
[~,sample_diagnostics]=clexpn(spin_system,@free_evolution,sample,'labframe');
spin_config=sample_diagnostics.ensemble.spin_config(1);
model=cce_compile(spin_system,sample,'labframe');
local=sample;
local.cce.bath_state=struct('mode','explicit_product',...
                            'rho',{spin_config.bath_rho});
local.cce.order=0;
small=clexpn(spin_system,@free_evolution,local,'labframe');
spin=pauli(2);
axes={'x','y','z'};
expectations=zeros(3,1);
for k=1:3
    expectations(k)=real(trace(spin.(axes{k})*spin_config.bath_rho{1}));
end
field=model.retained{1,2}*expectations;
x=[0 1;1 0]/2;
y=[0 -1i;1i 0]/2;
z=diag([1 -1])/2;
h=100*x+field(1)*x+field(2)*y+field(3)*z;
reference=zeros(numel(p.times),4);
for k=1:numel(p.times)
    u=expm(-1i*h*p.times(k));
    rho=u*p.rho0*u';
    reference(k,:)=rho(:).';
end
result=test_close(result,'static exterior exact',small,reference,1e-12,1e-12,...
                  'Independent observed-spin-only Hamiltonian with the same exterior assignment.');

% Exercise both explicitly named sampling estimators at exact full order
sample.output='density';
sample.cce.order=1;
sample.cce.averaging='after_expansion';
[after,after_diagnostics]=clexpn(spin_system,@free_evolution,sample,'labframe');
sample.cce.averaging='internal_then_expansion';
internal=clexpn(spin_system,@free_evolution,sample,'labframe');
result=test_close(result,'full-order sampling estimators',after,internal,...
                  1e-12,1e-12,'Exact raw density is linear in the initial state.');
result=test_true(result,'sample errors finite',all(isfinite(after_diagnostics.reduction.abs_error_bound),'all'),...
                 'Weighted sums retain explicit forward bounds.');

% Distinguish selective selected two-level transition pulses from full spin-one rotations

echo=sample;
echo.cce.bath_state=struct('mode','mixed_product');
echo.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
[observed,echo_diagnostics]=clexpn(spin_system,@hahn_echo,echo,'labframe');
model=cce_compile(spin_system,echo_diagnostics.parameters,'labframe');
cluster=cce_subsystem(model,2,struct(),[]);
[H,~,~,~,~]=cce_generators(cluster,echo_diagnostics.parameters,[]);
initial=kron(echo.rho0,eye(2)/2);
pulse=kron(-1i*[0 1;1 0],eye(2));
reference=zeros(size(observed));
for k=1:numel(echo.times)
    U=expm(-1i*H*echo.times(k)/2);
    if echo.times(k)>0
        U=U*pulse*U;
    end
    rho=U*initial*U';
    reduced=zeros(2);
    for a=1:2
        for b=1:2
            reduced(a,b)=trace(rho((a-1)*2+(1:2),(b-1)*2+(1:2)));
        end
    end
    reference(k,:)=reduced(:).';
end
result=test_close(result,'explicit pulse timestamps',observed,reference,1e-12,1e-12,...
                  'Independent propagation places the pulse at half the total elapsed time.');
a=[1;0;0];
b=[0;1;0];
sx=[0 1 0;1 0 1;0 1 0]/sqrt(2);
result=test_true(result,'spin-one selective distinction',...
                 norm(expm(-1i*pi*sx)-expm(-1i*pi*(a*b'+b*a')/2),'fro')>1,...
                 'A full physical spin-one pi rotation is not a selective transition pulse.');

metrics.trace_error=max(trace_error);
metrics.hermiticity_error=max(hermiticity_error);
metrics.minimum_eigenvalue=min(minimum_eigenvalue);
metrics.masked_elements=nnz(~ratio_reduced.valid);

% Report the measured metrics
disp(metrics);

end


