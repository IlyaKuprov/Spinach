% Test: kernel/gcce_observed_mediation | Observed-spin mediated gCCE
% Test bath-disconnected but observed-mediated gCCE correlations. Syntax:
%
%                  result=test_gcce_observed_mediation()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_observed_mediation()

% Demand the mandatory counterexample to bath-only connectivity pruning
result=new_test_result('kernel/gcce_observed_mediation','Observed-spin mediated gCCE',...
                       'Two bath spins with zero direct coupling share a flip-flopping centre.');
sys.isotopes={'1H','13C','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
spin_system=create(sys,inter);
spin_system=cce_test_write(spin_system,[1 2],diag([80 80 0]));
spin_system=cce_test_write(spin_system,[1 3],diag([110 110 0]));
p.observed=1;
p.cce.bath_spins=[2 3];
p.cce.approximation='generalized';
p.cce.expansion='standard';
p.output='signal';
p.normalized=true;
p.transition.spins=1;
p.transition.alpha=[1;0];
p.transition.beta=[0;1];
p.rho0=ones(2)/2;
p.orientation=[0 0 0];
p.cce.interaction_cutoff=Inf;
p.cce.graph=struct('mode','coupling','threshold',0,'score','coupling_norm');
p.cce.meanfield='none';
p.times=linspace(0,0.08,41);
p.cce.bath_state=struct('mode','mixed_product');
p.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Compare singletons with the retained exact observed-mediated pair
p.cce.order=1;
first=clexpn(spin_system,@free_evolution,p,'labframe');
p.cce.order=2;
[second,second_diagnostics]=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_true(result,'pair retained',numel(second_diagnostics.clusters)==3,...
                 'Bath-bath zero does not factorize transverse observed exchange.');

% Construct a separate full three-spin dense exponential reference
x=[0 1;1 0]/2;
y=[0 -1i;1i 0]/2;
H=80*(kron(x,kron(x,eye(2)))+kron(y,kron(y,eye(2))))+...
  110*(kron(x,kron(eye(2),x))+kron(y,kron(eye(2),y)));
rho=kron(ones(2)/2,eye(4)/4);
reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    u=expm(-1i*H*p.times(k));
    evolved=u*rho*u';
    reference(k)=2*trace(evolved(1:4,5:8));
end
result=test_close(result,'exact pair',second,reference,1e-11,1e-11,...
                  'Independent full three-spin Kronecker/expm reference.');
difference=max(abs(first-second));
result=test_true(result,'nonfactorization measured',difference>1e-2,...
                 'Observed-spin exchange creates a nonzero irreducible bath-pair factor.');
metrics.first_order_error=difference;
metrics.full_order_error=max(abs(second-reference));

% Resolve two independently coupled observed spins with joint detection
sys.isotopes={'1H','1H','13C','13C'};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);
spin_system=create(sys,inter);
spin_system=cce_test_write(spin_system,[1 3],diag([0 0 2]));
spin_system=cce_test_write(spin_system,[2 4],diag([0 0 2]));
p.observed=[1 2];
p=rmfield(p,'transition');
p.cce.bath_spins=[3 4];
p.rho0=ones(4)/4;
p.output='observable';
p.coil=kron(x,eye(2))+kron(eye(2),x);
p.times=[0 1];
p.cce.order=2;
reference=cos(p.times(:));
p.cce.averaging='internal_then_expansion';

% Check unitary and master-equation paths against the exact signal
for method={'gcce','megcce'}
    p.cce.approximation='generalized';
    p.cce.basis=struct('formalism','zeeman-hilb','approximation','none');
    if strcmp(method{1},'megcce')
        p.bath_jumps=struct([]);
        p.cce.basis.formalism='zeeman-liouv';
    end
    for mode={'mixed_product','sampled'}
        p.cce.bath_state=struct('mode',mode{1});
        if strcmp(mode{1},'sampled')
            p.cce.bath_state.samples=4;
            p.cce.bath_state.seed=11;
            p.cce.bath_state.basis='uniform';
            p.cce.bath_state.probabilities={ones(1,2)/2,ones(1,2)/2};
        end
        actual=clexpn(spin_system,@free_evolution,p,'labframe');
        name=['joint detection ' method{1} ' ' mode{1}];
        result=test_close(result,name,actual,reference,1e-8,1e-8,...
                          'Disconnected Hamiltonian components need not factorize a sum detection operator.');
        metrics.([method{1} '_' mode{1} '_joint_error'])=max(abs(actual-reference));
    end
    if isfield(p,'bath_jumps')
        p=rmfield(p,'bath_jumps');
    end
end

% Report the measured metrics
disp(metrics);

end


