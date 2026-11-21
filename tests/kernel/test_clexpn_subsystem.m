% Tests native compilation, retained ordering, configuration, and cache keys.
% Syntax:
%
%                   result=test_clexpn_subsystem()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_subsystem()

% State the physical domains and deterministic mixed-multiplicity fixture
result=new_test_result('kernel/clexpn_subsystem','Native CCE subsystems',...
                       'Canonical units, retained ordering, states and cache ownership.');
sys.isotopes={'E3','14N','13C'};
sys.magnet=0.035;
sys.output='hush';
inter.zeeman.matrix={diag([2.01 2.02 2.03]),diag([4 2 -3]),diag([1 -2 4])};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=[11 3 -2;4 17 5;7 -6 23];
inter.coupling.matrix{2,3}=[2 -1 3;4 5 1;-2 6 8];
inter.coupling.matrix{1,3}=diag([1 2 3]);
inter.coupling.matrix{2,2}=diag([-4 -4 8]);
inter.coordinates={[0 0 0],[1 0 0],[2 0 0]};
system=create(sys,inter);
parameters.observed=[3 1];
parameters.times=[0 1e-6];
parameters.orientation=[0.3 0.7 0.2];
parameters.cce.bath_spins=2;
parameters.cce.groups={2};
parameters.cce.backend='direct';
parameters.cce.numerics=struct('atol',1e-12,'rtol',1e-12,'max_memory_bytes',2e8);
model=cce_compile(system,parameters,'labframe');
[cluster,cache]=cce_subsystem(model,2,struct(),[]);
result=test_true(result,'Observed-first physical order',isequal(cluster.physical_ids,[3 1 2]),...
                 'Physical IDs must never be sorted before tensor embedding.');
result=test_true(result,'Mixed dimensions and units',...
                 isequal(cluster.multiplicities,[2 3 3])&&cluster.dimension==18&&...
                 isequal(cluster.coordinates_A,cell2mat(inter.coordinates([3 1 2]).'))&&...
                 isequal(cluster.isotopes,sys.isotopes([3 1 2]))&&...
                 isequal(cluster.groups,{2})&&~isfield(model.spin_system,'bas'),...
                 'Canonical isotope, geometry, group and dimension metadata must survive.');
result=test_true(result,'Bidirectional physical map',...
                 isequal(cluster.physical_to_local(cluster.local_to_physical),1:3),...
                 'Retained local indices must map back to the original physical IDs.');

% Reproduce the direct-versus-stock Hamiltonian fixture through native inputs
H=cce_hamiltonian(model,2,struct('backend','direct'));
[repeated,~]=cce_options(system,parameters,@free_evolution,'labframe');
reference=cce_compile(system,repeated,'labframe');
H_repeat=cce_hamiltonian(reference,2,struct('backend','direct'));
result=test_true(result,'Native compiler bitwise identity',isequal(H,H_repeat),...
                 'Repeated native compilation must preserve the retained numerical algebra.');
H_stock=cce_hamiltonian(model,2,struct('backend','stock'));
result=test_close(result,'Stock anisotropic reference',H,H_stock,5e-7,2e-12,...
                  'The canonical physical representation must agree with the independent reference.');

% Empty bath clusters retain every observed spin and the declared ordering
empty=cce_subsystem(model,[],struct(),[]);
result=test_true(result,'Empty bath retains observed space',...
                 isequal(empty.physical_ids,[3 1])&&empty.dimension==6&&...
                 empty.bath_dimension==1&&isempty(empty.bath_rho),...
                 'The empty cluster is the observed subsystem, not a scalar.');

% Cache signatures include physical model, assumptions, and representation
[again,cache]=cce_subsystem(model,2,struct(),cache);
result=test_true(result,'Exact cache hit',isequaln(cluster,again)&&cache.hits==1,...
                 'An unchanged request must reuse the same immutable view.');
changed=model;
changed.formalism='zeeman-liouv';
[again,cache]=cce_subsystem(changed,2,struct(),cache);
result=test_true(result,'Representation invalidates cache',cache.misses==2&&...
                 strcmp(again.formalism,'zeeman-liouv'),'Representation belongs in the signature.');
changed=cce_compile(system,parameters,'nmr');
[again,cache]=cce_subsystem(changed,2,struct(),cache);
result=test_true(result,'Assumptions invalidate cache',cache.misses==3&&...
                 isequal(again.assumptions,{'nmr'}),'Assumptions belong in the signature.');
changed.spin_system.inter.zeeman.matrix{1}(1,1)=1;
[~,cache]=cce_subsystem(changed,2,struct(),cache);
result=test_true(result,'Physical model invalidates cache',cache.misses==4,...
                 'Canonical model changes must never retain an old cache entry.');

% Nonuniform populations and exterior fields share one bath configuration
parameters.observed=1;
parameters.cce.bath_spins=[3 2];
parameters.cce.groups={[3 2]};
parameters.cce.meanfield='static';
parameters.cce.numerics=struct('atol',1e-12,'rtol',1e-12,'max_memory_bytes',2e8);
model=cce_compile(system,parameters,'labframe');
configuration.bath_ids=[2 3];
configuration.bath_rho={diag([0.2 0.3 0.5]),diag([0.8 0.2])};
[cluster,cache]=cce_subsystem(model,2,configuration,[]);
result=test_close(result,'Nonuniform local bath density',cluster.bath_rho{1},...
                  diag([0.2 0.3 0.5]),0,0,...
                  'The canonical physical representation must agree with the independent reference.');
result=test_close(result,'Exterior magnetic expectation',...
                  cluster.configuration.exterior_expectations,[0 0 0.3],1e-15,0,...
                  'The canonical physical representation must agree with the independent reference.');
configuration.bath_rho{2}=diag([0.1 0.9]);
[cluster,cache]=cce_subsystem(model,2,configuration,cache);
result=test_true(result,'Configuration invalidates cache',cache.misses==2,...
                 'Changing a physical bath state must invalidate the cached view.');
result=test_close(result,'Changed exterior field',...
                  cluster.configuration.exterior_expectations,[0 0 -0.4],1e-15,0,...
                  'The canonical physical representation must agree with the independent reference.');

% A scalar coupling has one Hz-to-angular-frequency conversion
clear sys inter parameters;
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
inter.zeeman.scalar={0,0};
inter.coupling.scalar={0,137;0,0};
system=create(sys,inter);
parameters.observed=1;
parameters.times=0;
parameters.cce.numerics=struct('atol',1e-12,'rtol',1e-12,'max_memory_bytes',2e8);
model=cce_compile(system,parameters,'labframe');
H=cce_hamiltonian(model,2,struct());
expected=2*pi*137*[1/4 0 0 0;0 -1/4 1/2 0;0 1/2 -1/4 0;0 0 0 1/4];
result=test_close(result,'Scalar Hz to rad/s',H,expected,1e-12,1e-12,...
                  'The canonical physical representation must agree with the independent reference.');

% Selected magnetic states are resolved before their projectors are compiled
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
model=cce_compile(system,parameters,'labframe');
result=test_close(result,'Resolved transition projector',...
                  model.manifold_projectors.alpha,diag([1 0]),0,0,...
                  'Native selectors must compile to the same ordered physical state.');
hilbert_hash=model.source_hash;

% Explicit representations remove stale whole-system basis arrays
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
system=basis(system,bas);
parameters.cce.numerics=struct('atol',1e-12,'rtol',1e-12,'max_memory_bytes',2e8);
model=cce_compile(system,parameters,'labframe');
result=test_true(result,'Carry representation without stale basis',...
                 strcmp(model.formalism,'zeeman-liouv')&&~isfield(model.spin_system,'bas')&&...
                 ~strcmp(model.source_hash,hilbert_hash),...
                 'Retained subsystems must not inherit full-system basis arrays.');

% Reject duplicate sites and invalid physical populations
rejected=false;
try
    cce_subsystem(model,[2 2],struct(),[]);
catch
    rejected=true;
end
result=test_true(result,'Duplicate physical site rejected',rejected,'Ownership must be unique.');
configuration=struct('bath_ids',2,'bath_rho',{{diag([1.1 -0.1])}});
rejected=false;
try
    cce_subsystem(model,2,configuration,[]);
catch
    rejected=true;
end
result=test_true(result,'Negative population rejected',rejected,'Physical populations must be positive.');

end


