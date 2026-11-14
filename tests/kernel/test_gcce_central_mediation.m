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
result.metrics.first_order_error=difference;
result.metrics.full_order_error=max(abs(second.value-reference));

end

