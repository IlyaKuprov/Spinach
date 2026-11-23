% Test: kernel/clexpn_descriptors | Native cluster descriptors
% Tests canonical cluster spin operators and explicit density descriptors.
% Syntax:
%
%                  result=test_clexpn_descriptors()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_descriptors()

% Mixed multiplicities and permuted observed and bath physical identifiers
result=new_test_result('kernel/clexpn_descriptors','Native cluster descriptors',...
                       'Stock normalisation, isotope sums, products and explicit populations.');
sys.isotopes={'1H','14N','13C','1H'};
sys.magnet=0;
sys.output='hush';
sys.tols.cce_atol=1e-12;
sys.tols.cce_rtol=1e-12;
sys.tols.cce_mem=2e8;
inter.zeeman.scalar={0,0,0,0};
system=create(sys,inter);
parameters.observed=[3 1];
parameters.times=0;
parameters.cce.bath_spins=[4 2];
model_p=cce_options(system,parameters,@cce_compile,'labframe');
model=cce_compile(system,model_p,'labframe');
[cluster,cache]=cce_subsystem(model,model_p,[4 2],struct(),[]);
[swapped,cache]=cce_subsystem(model,model_p,[2 4],struct(),cache);
result=test_true(result,'Retained ordering invalidates cache',cache.misses==2&&...
                 isequal(swapped.cce.ids,[3 1 2 4]),...
                 'Ordered physical IDs must be part of the exact cache signature.');
operators=pauli(3);
A=cce_operator(cluster,'Lz',2);
result=test_close(result,'Arbitrary-spin normalisation',A,...
                  kron(speye(8),operators.z),0,0,...
                  'The canonical physical representation must agree with the independent reference.');
A=cce_operator(cluster,{'Lx','Lz'},{3,2});
half=pauli(2);
result=test_close(result,'Observed-first tensor product',A,...
                  kron(kron(half.x,speye(4)),operators.z),0,0,...
                  'The canonical physical representation must agree with the independent reference.');
A=cce_operator(cluster,'Lz','1H');
expected=kron(kron(speye(2),half.z),speye(6))+...
         kron(kron(speye(4),half.z),speye(3));
result=test_close(result,'Isotope-group sum',A,expected,0,0,...
                  'The canonical physical representation must agree with the independent reference.');

% All supported descriptors agree with the stock Hilbert operator service
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
stock=basis(system,bas);
labels={'E','Lx','Ly','Lz','L+','L-'};
indices=reshape(1:24,fliplr([2 3 2 2]));
indices=permute(indices,[3 1 4 2]);
indices=indices(:);
for n=1:numel(labels)
    A=cce_operator(cluster,labels{n},2);
    reference=operator(stock,labels{n},2);
    reference=reference(indices,indices);
    result=test_close(result,['Stock descriptor ' labels{n}],A,reference,1e-14,1e-14,...
                  'The canonical physical representation must agree with the independent reference.');
end
rho=cce_state(cluster,{'Lx','Lz'},{3,2});
reference=state(stock,{'Lx','Lz'},{3,2});
reference=reference(indices,indices);
result=test_close(result,'Stock product state',rho,reference,1e-14,1e-14,...
                  'The canonical physical representation must agree with the independent reference.');

% Fully specified nonuniform populations are not normalised a second time
states={diag([0.7 0.3]),[0.2 0.3 0.5],[0.9 0.1],[0.4 0.6]};
rho=cce_state(cluster,states,[1 2 3 4]);
expected=kron(kron(kron(diag(states{3}),states{1}),diag(states{4})),diag(states{2}));
result=test_close(result,'Explicit physical state order',rho,expected,0,0,...
                  'The canonical physical representation must agree with the independent reference.');
result=test_close(result,'Explicit trace-one density',trace(rho),1,1e-14,0,...
                  'The canonical physical representation must agree with the independent reference.');

% Liouville operators and states use conventional column-major vectorisation
bas.formalism='zeeman-liouv';
stock=basis(system,bas);
model_p=cce_options(stock,parameters,@cce_compile,'labframe');
model=cce_compile(stock,model_p,'labframe');
cluster=cce_subsystem(model,model_p,[4 2],struct(),[]);
A=cce_operator(cluster,'Ly',2);
H=kron(speye(8),operators.y);
expected=kron(speye(24),H)-kron(H.',speye(24));
result=test_close(result,'Liouville commutator',A,expected,0,0,...
                  'The canonical physical representation must agree with the independent reference.');
vector=cce_state(cluster,states,[1 2 3 4]);
result=test_close(result,'Column-major physical density',vector,rho(:),0,0,...
                  'The canonical physical representation must agree with the independent reference.');
vector=cce_state(cluster,'Ly',2);
result=test_close(result,'Column-major operator state',vector,H(:),0,0,...
                  'The canonical physical representation must agree with the independent reference.');

% An empty bath still permits observed-spin preparation and detection
cluster=cce_subsystem(model,model_p,[],struct(),[]);
A=cce_operator(cluster,'Lx',3);
H=kron(half.x,speye(2));
result=test_close(result,'Empty-bath observed commutator',A,...
                  kron(speye(4),H)-kron(H.',speye(4)),0,0,...
                  'The canonical physical representation must agree with the independent reference.');

% Physical states cannot omit populations or refer to discarded spins
rejected=false;
try
    cce_state(cluster,{[0.4 0.6]},3);
catch
    rejected=true;
end
result=test_true(result,'Unprovided population rejected',rejected,...
                 'Every retained spin needs an explicit density for physical-state descriptors.');
rejected=false;
try
    cce_operator(cluster,'Lz',2);
catch
    rejected=true;
end
result=test_true(result,'Discarded physical ID rejected',rejected,...
                 'Local indices must never be mistaken for physical IDs.');
rejected=false;
try
    cce_operator(cluster,'Qz',3);
catch
    rejected=true;
end
result=test_true(result,'Unsupported descriptor rejected',rejected,...
                 'Unsupported spin-operator labels cannot silently change meaning.');

end


