% Test: kernel/cce_hamiltonian | CCE Hamiltonians
% Tests direct and stock CCE cluster Hamiltonians. Syntax:
%
%                    result=test_cce_hamiltonian()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_hamiltonian()

% State the tested physical domains
result=new_test_result('kernel/cce_hamiltonian','CCE Hamiltonians',...
                       'Check arbitrary spin, ordering, rotations, retention and fallback.');

% Create a deterministic mixed-spin anisotropic fixture
sys.isotopes={'E3','14N','13C'};
sys.magnet=0.035;
sys.output='hush';
inter.zeeman.matrix={diag([2.01 2.02 2.03]),diag([4 2 -3]),diag([1 -2 4])};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=[11 3 -2;4 17 8;6 -5 23]*100;
inter.coupling.matrix{2,3}=[8 -5 4;3 -7 9;1 6 13]*31;
inter.coupling.matrix{1,3}=diag([7 9 13])*17;
inter.coupling.matrix{2,2}=[2 1 0;1 -3 2;0 2 1]*41;
system=create(sys,inter);
parameters.times=0;
parameters.observed=2;
parameters.cce.bath_spins=[3 1];
parameters.orientation=[0 0 0];
parameters.cce.interaction_cutoff=Inf;
parameters.cce.numerics=struct('atol',1e-5,'rtol',1e-12,'max_memory_bytes',1e8);

% Exercise every validated assumption, Euler choice, and retention selection
max_error=0;
max_relative=0;
max_propagator=0;
comparisons=0;
for angles={[0 0 0],[0.3 0.7 -0.2]}
    parameters.orientation=angles{1};
    for assumption={'labframe','nmr','esr','deer','deer-zz'}
        for retention={'','zeeman','couplings'}
            assumptions=assumption{1};
            if ~isempty(retention{1})
                assumptions={assumption{1},retention{1}};
            end
            model=cce_compile(system,parameters,assumptions);
            [H,meta]=cce_hamiltonian(model,[3 1],struct('backend','direct'));
            reference=cce_hamiltonian(model,[3 1],struct('backend','stock'));
            error_norm=norm(H-reference,'fro');
            max_error=max(max_error,error_norm);
            max_relative=max(max_relative,error_norm/max(norm(reference,'fro'),realmin));
            comparisons=comparisons+1;
            propagator_error=norm(expm(full(-1i*H*1e-7))-...
                                  expm(full(-1i*reference*1e-7)),'fro');
            max_propagator=max(max_propagator,propagator_error);
            result=test_true(result,'dense propagator parity',...
                             propagator_error<=parameters.cce.numerics.atol,...
                             'Independent dense exponentials obey the requested observable budget.');
            result=test_true(result,[assumption{1} '/' retention{1} ' parity'],...
                             error_norm<=meta.roundoff_bound&&error_norm<=meta.observable_budget,'Stock/direct error must satisfy roundoff and requested budgets.');
        end
    end
end
fprintf('CCE Hamiltonian: %d parity cases, max Frobenius=%g, max relative=%g\n',...
        comparisons,max_error,max_relative);
fprintf('CCE Hamiltonian: max dense propagator error=%g\n',max_propagator);

% Demonstrate that transposing a single ordered nonsymmetric cell is detected
assumptions='labframe';
parameters.orientation=[0 0 0];
model=cce_compile(system,parameters,assumptions);
[reference,meta]=cce_hamiltonian(model,[3 1],struct('backend','stock'));
mutant=system;
mutant.inter.coupling.matrix{1,2}=mutant.inter.coupling.matrix{1,2}.';
wrong=cce_hamiltonian(cce_compile(mutant,parameters,assumptions),[3 1],struct());
result=test_true(result,'transpose mutation rejected',...
                 norm(wrong-reference,'fro')>100*meta.roundoff_bound,'Antisymmetric ordered-cell convention is observable.');

% Check one-sided versus two-half-cell storage for a nonsymmetric tensor
half=system;
tensor=half.inter.coupling.matrix{1,2};
half.inter.coupling.matrix{1,2}=tensor/2;
half.inter.coupling.matrix{2,1}=tensor.'/2;
H=cce_hamiltonian(cce_compile(half,parameters,assumptions),[3 1],struct());
result=test_close(result,'nonsymmetric half cells',...
                  H,reference,meta.roundoff_bound,0,'Reverse ordered cell requires tensor transpose.');

% Check covariance under a physical relabelling independent of product order
permutation=[3 1 2];
permuted_sys=sys;
permuted_sys.isotopes=sys.isotopes(permutation);
permuted_inter=inter;
permuted_inter.zeeman.matrix=inter.zeeman.matrix(permutation);
permuted_inter.coupling.matrix=inter.coupling.matrix(permutation,permutation);
permuted_params=parameters;
permuted_params.observed=3;
permuted_params.cce.bath_spins=[1 2];
permuted=cce_compile(create(permuted_sys,permuted_inter),permuted_params,assumptions);
H=cce_hamiltonian(permuted,[1 2],struct());
result=test_close(result,'physical spin permutation',...
                  H,reference,meta.roundoff_bound,0,'Relabelling leaves the ordered physical Hamiltonian invariant.');

% Check an observed-spin system in nonascending physical order and a proper subset
parameters.observed=[3 1];
parameters.cce.bath_spins=2;
model=cce_compile(system,parameters,assumptions);
[H,meta]=cce_hamiltonian(model,2,struct());
reference=cce_hamiltonian(model,2,struct('backend','stock'));
result=test_close(result,'multiple observed spins',...
                  H,reference,meta.roundoff_bound,0,'Product basis follows the supplied observed order.');
[H,meta]=cce_hamiltonian(model,[],struct());
reference=cce_hamiltonian(model,[],struct('backend','stock'));
result=test_close(result,'observed-spin-only subset',...
                  H,reference,meta.roundoff_bound,0,'Metadata slicing preserves physical indices.');

% Compare exterior fields with an explicit partial expectation of the full system
parameters.observed=1;
parameters.cce.bath_spins=[2 3];
model=cce_compile(system,parameters,assumptions);
[full_h,full_meta]=cce_hamiltonian(model,[2 3],struct());
spin_config.exterior_ids=3;
spin_config.exterior_expectations=[0 0 1/2];
[H,meta]=cce_hamiltonian(model,2,spin_config);
embedding=kron(speye(meta.dimension),[1;0]);
reference=embedding'*full_h*embedding;
reference=reference-trace(reference-H)*speye(meta.dimension)/meta.dimension;
result=test_close(result,'exterior expectation',...
                  H,reference,full_meta.roundoff_bound,0,'Exterior self energy is a removable scalar.');

% Test finite physical cutoffs independently of graph pruning
cut=system;
cut.inter.coordinates={[0 0 0],[0 0 2],[0 0 20]};
parameters.cce.interaction_cutoff=3;
model=cce_compile(cut,parameters,assumptions);
[H,meta]=cce_hamiltonian(model,[2 3],struct());
reference=cce_hamiltonian(model,[2 3],struct('backend','stock'));
result=test_close(result,'physical distance cutoff',...
                  H,reference,meta.roundoff_bound,0,'Both paths omit the same declared pairs.');
result=test_true(result,'cutoff omission count',...
                 size(model.omitted_pairs,1)==2,'Two physical pairs exceed three Angstrom.');

% Reject a declared memory budget before allocating cluster-sized matrices
small=parameters;
small.cce.interaction_cutoff=Inf;
small.cce.numerics.max_memory_bytes=1;
rejected=false;
try
    cce_hamiltonian(cce_compile(system,small,assumptions),[2 3],struct());
catch exception
    rejected=contains(exception.message,'max_memory_bytes');
end
result=test_true(result,'memory ceiling rejected',...
                 rejected,'Insufficient explicit memory must stop before operator allocation.');

% Reject an assumption name that has no declared physical meaning
parameters.cce.interaction_cutoff=Inf;
assumptions='unsupported';
rejected=false;
try
    cce_compile(system,parameters,assumptions);
catch
    rejected=true;
end
result=test_true(result,'unsupported assumption rejected',...
                 rejected,'A physical approximation cannot be inferred from an unknown name.');
assumptions='labframe';

% Force automatic fallback for a nonzero rank-three giant-spin term
sys.isotopes={'E4'};
sys.magnet=0.01;
inter=struct();
inter.giant.coeff={{zeros(3,1),zeros(5,1),[0;0;0;37;0;0;0]}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0]}};
parameters.observed=1;
parameters.cce.bath_spins=[];
parameters.cce.interaction_cutoff=Inf;
model=cce_compile(create(sys,inter),parameters,assumptions);
[H,meta]=cce_hamiltonian(model,[],struct());
reference=cce_hamiltonian(model,[],struct('backend','stock'));
result=test_close(result,'rank-three fallback',...
                  H,reference,meta.roundoff_bound,0,'Automatic fallback retains giant-spin terms.');
result=test_true(result,'fallback recorded',...
                 strcmp(meta.backend,'stock'),'Unsupported direct domains cannot be guessed.');
rejected=false;
try
    cce_hamiltonian(model,[],struct('backend','direct'));
catch
    rejected=true;
end
result=test_true(result,'unsupported direct rejection',...
                 rejected,'Explicit direct must fail closed.');

end


