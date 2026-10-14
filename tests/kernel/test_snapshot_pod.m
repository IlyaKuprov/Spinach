% Tests weighted snapshot POD, complex affine Galerkin projection, and
% reduced observables against known subspaces, a triangular dissipative
% solution, and ordinary Spinach evolution. Syntax:
%
%                      result=test_snapshot_pod()
%
% Outputs:
%
%     result - maintained regression record with independent references
%
% Compact well-scaled identities use 1e-12 absolute tolerances. Snapshot
% residuals diagnose the supplied data, not the observable approximation.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il

function result=test_snapshot_pod()

% Initialise the maintained test record and preserve the caller's random state
result=new_test_result('kernel/snapshot_pod','Snapshot POD model reduction',...
                      'Explicit weighted subspaces and projected observables must reproduce independent references.');
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng));
exact_opts=struct('method','svd','oversampling',0,'power_iters',0,'seed',7);
rand_opts=struct('method','randomised','oversampling',1,'power_iters',2,'seed',7);

% Test a weighted complex diagonal snapshot matrix with known singular values
blocks={diag([3 1i 0.5]);diag([1 2i 1])}; weights=[1;4];
observed=pod_basis(blocks,weights,2,exact_opts);
result=test_close(result,'weighted singular values',observed.singular_vals,...
                  sqrt([17;13;4.25]),1e-12,0,'block-energy weights enter through their square roots');
result=test_close(result,'weighted retained projector',observed.basis*observed.basis',...
                  diag([1 1 0]),1e-12,0,'the two largest weighted modes span the first two axes');
result=test_close(result,'absolute residual',observed.residual_abs,sqrt(4.25),...
                  1e-12,0,'the omitted third mode supplies the Frobenius residual');
result=test_close(result,'relative residual',observed.residual_rel,sqrt(4.25/34.25),...
                  1e-12,0,'relative residual divides by the supplied weighted data norm');

% Check factors, zero weights, and the full-rank projection diagnostic
raw=[1 2i 3;2 -1i 1;1i 0 2]; [left,values,~]=svd(raw,'econ');
raw_basis=pod_basis({raw},1,2,exact_opts);
factor_basis=pod_basis({left*values},1,2,exact_opts);
result=test_close(result,'equivalent factors',raw_basis.basis*raw_basis.basis',...
                  factor_basis.basis*factor_basis.basis',1e-12,0,...
                  'U*S factors with the same covariance give the same POD subspace');
zero_weight=pod_basis({raw;ones(3)},[1;0],3,exact_opts);
result=test_close(result,'zero-weight ignored block',zero_weight.residual_abs,0,...
                  1e-12,0,'a full-rank basis represents the nonzero-weight data exactly');

% Check seeded sketches and random-state restoration on success and rank failure
rng(123,'twister'); before=rng;
random_a=pod_basis({raw},1,2,rand_opts); after=rng;
random_b=pod_basis({raw},1,2,rand_opts);
result=test_true(result,'caller RNG preserved',isequal(before,after)&&isequal(before,rng),...
                 'seeded construction must not change subsequent caller random draws');
result=test_close(result,'seed repeatability',random_a.basis,random_b.basis,0,0,...
                  'the same explicit seed and inputs reproduce the sketch');
result=test_close(result,'full sketch subspace',random_a.basis*random_a.basis',...
                  raw_basis.basis*raw_basis.basis',1e-12,0,...
                  'a sketch spanning the full data space recovers the exact retained subspace');
failed=false;
try
    pod_basis({diag([1 0 0])},1,2,rand_opts);
catch
    failed=true;
end
result=test_true(result,'rank failure restores RNG',failed&&isequal(before,rng),...
                 'rank-deficient requests fail without supplying null-space directions or changing RNG');

% Compare genuinely truncated blockwise sketches with a concatenated reference
rng(632,'twister'); snapshots=complex(randn(18,27),randn(18,27));
split_blocks={sparse(snapshots(:,1:7));snapshots(:,8:19);sparse(snapshots(:,20:end))};
split_weights=[0.5;2;0];
weighted=[sqrt(0.5)*snapshots(:,1:7) sqrt(2)*snapshots(:,8:19) zeros(18,8)];
for power_iters=[0 2]
    sketch_opts=struct('method','randomised','oversampling',2,...
                       'power_iters',power_iters,'seed',17);
    rng(17,'twister'); omega=complex(randn(27,6),randn(27,6));
    [range,~]=qr(weighted*omega,0);
    for k=1:power_iters
        [dual_range,~]=qr(weighted'*range,0);
        [range,~]=qr(weighted*dual_range,0);
    end
    [left_modes,values]=svd(range'*weighted,'econ');
    reference=range*left_modes(:,1:4); before=rng;
    observed=pod_basis(split_blocks,split_weights,4,sketch_opts);
    result=test_close(result,['blockwise projector ' num2str(power_iters)],...
                      observed.basis*observed.basis',reference*reference',1e-12,0,...
                      'mixed sparse and dense weighted blocks recover the explicit truncated sketch');
    result=test_close(result,['blockwise values ' num2str(power_iters)],...
                      observed.singular_vals,diag(values),1e-12,0,...
                      'compressed singular values match the independently concatenated reference');
    result=test_close(result,['blockwise residual ' num2str(power_iters)],...
                      observed.residual_abs,norm(weighted-observed.basis*(observed.basis'*weighted),'fro'),...
                      1e-12,0,'the blockwise residual retains all supplied weights');
    result=test_true(result,['blockwise RNG ' num2str(power_iters)],isequal(before,rng),...
                     'genuinely truncated sketches preserve the caller random state');
end

% Check sparse snapshots against dense data and reject nonfinite or zero data
for sparse_opts={exact_opts,rand_opts}
    sparse_result=pod_basis({sparse(raw);sparse(3,2)},[1;0],2,sparse_opts{1});
    dense_result=pod_basis({raw;zeros(3,2)},[1;0],2,sparse_opts{1});
    result=test_close(result,['sparse snapshots ' sparse_opts{1}.method],...
                      sparse_result.basis*sparse_result.basis',...
                      dense_result.basis*dense_result.basis',1e-12,0,...
                      'finite sparse snapshots and zero-weight blocks preserve the weighted subspace');
    for bad_value=[NaN Inf 0]
        before=rng; rejected=false;
        try
            pod_basis({sparse(1,1,bad_value,3,3)},1,1,sparse_opts{1});
        catch
            rejected=true;
        end
        result=test_true(result,['sparse data ' sparse_opts{1}.method ' ' num2str(bad_value)],...
                         rejected&&isequal(before,rng),...
                         'stored NaN/Inf and all-zero sparse data fail without changing caller RNG');
    end
end

% Retain sparse Liouvillian terms and reject nonfinite stored entries
sparse_basis=[speye(2);sparse(2046,2)]; sparse_state=sparse_basis(:,1);
sparse_model=pod_project(full(sparse_basis),{speye(2048)},full(sparse_state),full(sparse_state));
result=test_close(result,'sparse term projection',sparse_model.terms{1},eye(2),...
                  1e-12,0,'sparse implicit zeros do not alter the Galerkin projection');
zero_model=pod_project(full(sparse_basis),{sparse(2048,2048)},...
                       full(sparse_state),full(sparse_state));
result=test_close(result,'zero sparse term',zero_model.terms{1},zeros(2),...
                  0,0,'a finite zero Liouvillian is a valid affine term');
for bad_value=[NaN Inf]
    rejected=false;
    try
        pod_project(full(sparse_basis),{sparse(3,3,bad_value,2048,2048)},...
                    full(sparse_state),full(sparse_state));
    catch
        rejected=true;
    end
    result=test_true(result,['sparse nonfinite ' num2str(bad_value)],rejected,...
                     'sparse validation still rejects nonfinite stored values');
end

% Solve a non-Hermitian triangular dissipative model in an invariant subspace
H=diag([2 -3 7]); R=[-0.4 0.2 0;0 -0.8 0;0 0 -1.2];
rho=[1+1i;-0.3+0.2i;0]; coil=[0.4-0.2i;1+1i;2];
subspace=[1 1i;1i 1;0 0]/sqrt(2); terms={H;1i*R};
model=pod_project(subspace,terms,rho,coil); coefficients=[0.7;1.3];
L=coefficients(1)*H+1i*coefficients(2)*R;
result=test_close(result,'affine projection',...
                  coefficients(1)*model.terms{1}+coefficients(2)*model.terms{2},...
                  subspace'*L*subspace,1e-12,0,'fixed-basis projection commutes with supplied affine assembly');
times=(0:20)*0.03; rates=-1i*diag(L);
second=rho(2)*exp(rates(2)*times);
first=rho(1)*exp(rates(1)*times)+coefficients(2)*R(1,2)*rho(2)*...
      (exp(rates(2)*times)-exp(rates(1)*times))/(rates(2)-rates(1));
reference=conj(coil(1))*first+conj(coil(2))*second;
fid=pod_evolve(model,coefficients,0.03,20);
result=test_close(result,'invariant dissipative trace',fid,reference,1e-12,0,...
                  'the independent triangular solution checks damping, the -1i sign, and detection conjugation');
full_model=pod_project(eye(3),terms,rho,coil);
result=test_close(result,'full-rank dissipative trace',...
                  pod_evolve(full_model,coefficients,0.03,20),reference,1e-12,0,...
                  'full-rank Galerkin dynamics is a coordinate change only');
result=test_close(result,'initial and zero-time trace',pod_evolve(model,coefficients,0,2),...
                  repmat(coil'*rho,1,3),1e-12,0,'the initial sample and zero-time repetitions retain the complex observable');

% Compare a damped proton coherence and longitudinal component with Spinach
sys.magnet=14.1; sys.isotopes={'1H'};
inter.zeeman.scalar={0}; inter.relaxation={'t1_t2'};
inter.r1_rates={0.7}; inter.r2_rates={1.1}; inter.equilibrium='zero'; inter.rlx_keep='labframe';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
spin_system.tols.prop_chop=eps; spin_system.tols.liouv_zero=eps;
spin_system.tols.zte_tol=eps;
H=2*pi*0.8*operator(spin_system,'Lz',1); R=relaxation(spin_system);
rho=state(spin_system,'L+',1)+0.3i*state(spin_system,'Lz',1); coil=rho;
trajectory=evolution(spin_system,H+1i*R,[],rho,0.04,20,'trajectory');
phys_basis=pod_basis({full(trajectory)},1,2,exact_opts);
phys_model=pod_project(phys_basis.basis,{H;1i*R},full(rho),full(coil));
phys_fid=pod_evolve(phys_model,[1;1],0.04,20);
spinach_fid=evolution(spin_system,H+1i*R,coil,rho,0.04,20,'observable');
result=test_close(result,'actual Spinach evolution',phys_fid,spinach_fid.',1e-12,0,...
                  'complex Liouville states and dissipative relaxation match the maintained full evolution');

% Exercise public input contracts rather than silently adapting shapes or ranks
bad_model=model; bad_model.rho=bad_model.rho.';
invalid={@()pod_basis(blocks,[0;0],1,exact_opts),...
         @()pod_basis({zeros(3)},1,1,exact_opts),...
         @()pod_basis({eye(3);eye(2)},[1;1],1,exact_opts),...
         @()pod_basis({raw},-1,1,exact_opts),...
         @()pod_basis({raw},1,4,exact_opts),...
         @()pod_project(2*subspace,terms,rho,coil),...
         @()pod_project(subspace,{eye(2)},[1;0;0],[1;0;0]),...
         @()pod_evolve(model,[1 1],0.1,2),...
         @()pod_evolve(bad_model,coefficients,0.1,2),...
         @()pod_evolve(model,coefficients,-0.1,2),...
         @()pod_evolve(model,coefficients,0.1,1.5)};
for n=1:numel(invalid)
    rejected=false;
    try
        invalid{n}();
    catch
        rejected=true;
    end
    result=test_true(result,['invalid input ' num2str(n)],rejected,...
                     'invalid data, coordinates, and time contracts must be rejected');
end
clear restore_rng

end

