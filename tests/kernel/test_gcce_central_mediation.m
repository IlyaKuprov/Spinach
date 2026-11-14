% Test bath-disconnected but central-mediated gCCE correlations. Syntax:
%
%                  result=test_gcce_central_mediation()
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_central_mediation()

% Demand the mandatory counterexample to bath-only connectivity pruning
result=new_test_result('kernel/gcce_central_mediation','Central mediated gCCE',...
                       'Two bath spins with zero direct coupling share a flip-flopping centre.');
sys.isotopes={'1H','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)}; inter.coupling.matrix=cell(3);
spin_system=create(sys,inter);
spin_system.inter.coupling.matrix{1,2}=diag([80 80 0]);
spin_system.inter.coupling.matrix{1,3}=diag([110 110 0]);
p.central_spins=1; p.bath_spins=[2 3]; p.method='gcce'; p.expansion='standard';
p.quantity='coherence'; p.normalized=true; p.central.alpha=[1;0]; p.central.beta=[0;1];
p.central.rho=ones(2)/2; p.assumptions='labframe'; p.orientation=[0 0 0];
p.interaction_cutoff=Inf; p.graph=struct('mode','coupling','threshold',0,'score','coupling_norm');
p.meanfield='none'; p.times=linspace(0,0.08,41); p.sequence=[];
p.bath_state=struct('mode','mixed_product');
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Compare singletons with the retained exact central-mediated pair
p.order=1; first=cce(spin_system,p); p.order=2; second=cce(spin_system,p);
result=test_true(result,'pair retained',numel(second.clusters)==3,...
                 'Bath-bath zero does not factorize transverse central exchange.');

% Construct a separate full three-spin dense exponential reference
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2;
H=80*(kron(x,kron(x,eye(2)))+kron(y,kron(y,eye(2))))+...
  110*(kron(x,kron(eye(2),x))+kron(y,kron(eye(2),y)));
rho=kron(ones(2)/2,eye(4)/4); reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    u=expm(-1i*H*p.times(k)); evolved=u*rho*u';
    reference(k)=2*trace(evolved(1:4,5:8));
end
result=test_close(result,'exact pair',second.value,reference,1e-11,1e-11,...
                  'Independent full three-spin Kronecker/expm reference.');
difference=max(abs(first.value-second.value));
result=test_true(result,'nonfactorization measured',difference>1e-2,...
                 'Central exchange creates a nonzero irreducible bath-pair factor.');
metrics.first_order_error=difference;
metrics.full_order_error=max(abs(second.value-reference));

% Resolve two independently coupled central spins with a joint readout
sys.isotopes={'1H','1H','13C','13C'};
inter.zeeman.matrix=repmat({zeros(3)},1,4); inter.coupling.matrix=cell(4);
spin_system=create(sys,inter);
spin_system.inter.coupling.matrix{1,3}=diag([0 0 2]);
spin_system.inter.coupling.matrix{2,4}=diag([0 0 2]);
p.central_spins=[1 2]; p.bath_spins=[3 4]; p.central=struct('rho',ones(4)/4);
p.quantity='observable'; p.operators=kron(x,eye(2))+kron(eye(2),x);
p.times=[0 1]; p.order=2; reference=cos(p.times(:));
p.averaging='internal_then_expansion';

% Check unitary and master-equation paths against the exact signal
for method={'gcce','megcce'}
    p.method=method{1};
    if strcmp(p.method,'megcce'), p.jumps=struct([]); end
    for mode={'mixed_product','sampled'}
        p.bath_state=struct('mode',mode{1});
        if strcmp(mode{1},'sampled')
            p.bath_state.samples=4; p.bath_state.seed=11; p.bath_state.basis='uniform';
            p.bath_state.probabilities={ones(2,1)/2,ones(2,1)/2};
        end
        actual=cce(spin_system,p);
        name=['joint readout ' method{1} ' ' mode{1}];
        result=test_close(result,name,actual.value,reference,1e-8,1e-8,...
                          'Disconnected Hamiltonian components need not factorize a sum readout.');
        metrics.([method{1} '_' mode{1} '_joint_error'])=max(abs(actual.value-reference));
    end
    if isfield(p,'jumps'), p=rmfield(p,'jumps'); end
end

% Report the measured metrics
disp(metrics);

end


