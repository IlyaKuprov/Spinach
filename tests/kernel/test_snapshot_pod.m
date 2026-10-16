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

% Check tiny nonzero norms without using underflowed energy as a zero test
tiny_opts={exact_opts,...
           struct('method','randomised','oversampling',0,'power_iters',0,'seed',7),...
           struct('method','randomised','oversampling',0,'power_iters',2,'seed',7)};
for n=1:numel(tiny_opts)
    opts=tiny_opts{n}; before=rng; label=num2str(n); scale=1e-200;
    tiny=pod_basis({scale},1,1,opts);
    result=test_close(result,['tiny scalar value ' label],tiny.singular_vals/scale,...
                      1,1e-12,0,'the nonzero rank-one singular value survives energy underflow');
    result=test_close(result,['tiny scalar residual ' label],tiny.residual_rel,...
                      0,1e-12,0,'relative residual uses the unsquared norm');
    result=test_true(result,['tiny energy and RNG ' label],...
                     tiny.total_energy==0&&isequal(before,rng),...
                     'rounded-zero energy is a diagnostic and the caller RNG is preserved');
    if strcmp(opts.method,'randomised'), opts.oversampling=1; end
    tiny=pod_basis({scale*diag([3 1 0.5])},1,2,opts);
    result=test_close(result,['tiny truncated projector ' label],...
                      tiny.basis*tiny.basis',diag([1 1 0]),1e-12,0,...
                      'tiny rank-two truncation retains the independently known subspace');
    result=test_close(result,['tiny truncated residual ' label],tiny.residual_abs/scale,...
                      0.5,1e-12,0,'scaled absolute residual retains the omitted mode');
    result=test_close(result,['tiny relative residual ' label],tiny.residual_rel,...
                      0.5/sqrt(10.25),1e-12,0,'relative diagnostics remain meaningful at tiny scale');
    invalid_tiny={zeros(3),scale*diag([1 0 0])};
    for k=1:numel(invalid_tiny)
        rejected=false; before=rng;
        try
            pod_basis(invalid_tiny(k),1,2,opts);
        catch
            rejected=true;
        end
        result=test_true(result,['tiny invalid ' label ' ' num2str(k)],...
                         rejected&&isequal(before,rng),...
                         'zero data and unsupported rank fail with RNG preserved');
    end
end

% Check scale covariance without using physical squared energy for rank
for method={'svd','randomised'}
    opts=exact_opts; opts.method=method{1};
    if strcmp(opts.method,'randomised'), opts.oversampling=1; end
    for magnitude=[1e-200 1 1e200]
        data=magnitude*diag([3 1 0.5]); before=rng;
        observed=pod_basis({data},1,2,opts);
        result=test_close(result,['scale singular values ' method{1} ' ' num2str(magnitude)],...
                          observed.singular_vals/magnitude,[3;1;0.5],1e-12,0,...
                          'binary-normalised factorisation retains physical singular values');
        result=test_close(result,['scale residual ' method{1} ' ' num2str(magnitude)],...
                          observed.residual_rel,0.5/sqrt(10.25),1e-12,0,...
                          'relative reconstruction error is independent of physical amplitude');
        result=test_true(result,['scale RNG ' method{1} ' ' num2str(magnitude)],isequal(before,rng),...
                         'scale-normalised exact and seeded randomised paths preserve caller RNG');
    end
end

% Combine compensating weights and block amplitudes without large intermediates
for method={'svd','randomised'}
    opts=exact_opts; opts.method=method{1};
    if strcmp(opts.method,'randomised'), opts.oversampling=1; end
    data={2^-511*diag([3 1]);2^511*diag([1 2])};
    observed=pod_basis(data,[2^1022;2^-1022],1,opts);
    result=test_close(result,['compensated weights ' method{1}],observed.singular_vals,...
                      sqrt([10;5]),1e-12,0,'opposing block and weight scales define the same weighted operator');
    result=test_close(result,['compensated residual ' method{1}],observed.residual_rel,...
                      sqrt(1/3),1e-12,0,'the independent weighted Frobenius residual is retained');
    opts.oversampling=0; observed=pod_basis({realmax},realmax,1,opts);
    [fraction,~]=log2(sqrt(realmax)); [block_fraction,~]=log2(realmax);
    result=test_close(result,['scaled diagnostic frame ' method{1}],observed.singular_vals_scaled,...
                      fraction*block_fraction,1e-12,0,'the finite binary diagnostic frame survives physical overflow');
    result=test_true(result,['physical diagnostic overflow ' method{1}],...
                     isinf(observed.singular_vals)&&isinf(observed.total_energy)&&observed.residual_rel<1e-12,...
                     'physical diagnostic overflow does not invalidate a well-defined POD subspace');
end

% Reject unsupported training dynamic range without changing the caller RNG
before=rng; rejected=false;
try
    pod_basis({diag([realmax 2^-1074])},1,1,exact_opts);
catch err
    rejected=contains(err.message,'dynamic range');
end
result=test_true(result,'training range rejection',rejected&&isequal(before,rng),...
                 'normalisation never silently deletes a nonzero supplied component');

% Check that the Gram budget accepts native complex producer outputs
accepted=true;
for n=1:100
    sample=randn+1i*randn; produced=pod_basis({sample},1,1,exact_opts);
    try
        pod_project(produced.basis,{1},1,1);
    catch
        accepted=false;
    end
end
result=test_true(result,'native complex producer Gram budget',accepted,...
                 'dimension-dependent complex roundoff admits the native scalar SVD basis');

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
weighted=[sqrt(0.5)*snapshots(:,1:7) sqrt(2)*snapshots(:,8:19)];
for power_iters=[0 2]
    sketch_opts=struct('method','randomised','oversampling',2,...
                       'power_iters',power_iters,'seed',17);
    rng(17,'twister'); omega=complex(randn(19,6),randn(19,6));
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

% Check large oversampling without changing the retained subspace or spectrum
wide_data=[diag(32:-1:1) zeros(32)];
wide_opts=struct('method','randomised','oversampling',30,'power_iters',0,'seed',17);
wide_basis=pod_basis({sparse(wide_data)},1,2,wide_opts);
result=test_close(result,'oversampled retained projector',wide_basis.basis*wide_basis.basis',...
                  diag([1 1 zeros(1,30)]),1e-12,0,...
                  'a full-row sketch retains only the two leading diagonal snapshot directions');
result=test_close(result,'oversampled singular values',wide_basis.singular_vals,(32:-1:1)',...
                  1e-12,0,'retained mode construction does not discard compressed singular values');

% Compare real-only sketches with a real Gaussian reference and retain real storage
rng(632,'twister'); real_snapshots=randn(18,27);
real_blocks={sparse(real_snapshots(:,1:7));real_snapshots(:,8:19);1i*ones(18,8)};
real_weights=[0.5;2;0]; weighted=[sqrt(0.5)*real_snapshots(:,1:7) sqrt(2)*real_snapshots(:,8:19)];
for power_iters=[0 2]
    sketch_opts=struct('method','randomised','oversampling',2,...
                       'power_iters',power_iters,'seed',17);
    rng(17,'twister'); omega=randn(19,6); [range,~]=qr(weighted*omega,0);
    for k=1:power_iters
        [dual_range,~]=qr(weighted'*range,0);
        [range,~]=qr(weighted*dual_range,0);
    end
    [left_modes,values]=svd(range'*weighted,'econ');
    reference=range*left_modes(:,1:4); before=rng;
    observed=pod_basis(real_blocks,real_weights,4,sketch_opts);
    repeated=pod_basis(real_blocks,real_weights,4,sketch_opts);
    result=test_close(result,['real projector ' num2str(power_iters)],...
                      observed.basis*observed.basis',reference*reference',1e-12,0,...
                      'real weighted blocks recover the independent real truncated sketch');
    result=test_close(result,['real values ' num2str(power_iters)],observed.singular_vals,...
                      diag(values),1e-12,0,'real compressed singular values match the concatenated reference');
    result=test_close(result,['real residual ' num2str(power_iters)],observed.residual_abs,...
                      norm(weighted-observed.basis*(observed.basis'*weighted),'fro'),1e-12,0,...
                      'a zero-weight complex block does not contribute to the real projection residual');
    result=test_true(result,['real storage ' num2str(power_iters)],isreal(observed.basis),...
                     'complex storage is unnecessary when every positively weighted block is real');
    result=test_true(result,['real RNG ' num2str(power_iters)],isequal(before,rng),...
                     'real truncated sketches preserve the caller random state');
    result=test_close(result,['real repeatability ' num2str(power_iters)],...
                      observed.basis,repeated.basis,0,0,'real sketches are reproducible for a fixed seed');
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

% Preserve spectral padding when zero-weight columns exceed the active dimensions
thin_data=[3 0;0 1;0 0;0 0]; excluded=sparse(4,10000);
for method={'svd','randomised'}
    if strcmp(method{1},'svd')
        thin_opts=exact_opts;
    else
        thin_opts=struct('method','randomised','oversampling',2,'power_iters',2,'seed',7);
    end
    observed=pod_basis({thin_data;excluded},[1;0],2,thin_opts);
    result=test_close(result,['excluded spectrum ' method{1}],observed.singular_vals,...
                      [3;1;0;0],1e-12,0,'omitted zero-weight directions remain explicit zero singular values');
    result=test_close(result,['excluded subspace ' method{1}],observed.basis*observed.basis',...
                      diag([1 1 0 0]),1e-12,0,'oversampling beyond the active columns cannot add data directions');
    if strcmp(method{1},'randomised'), thin_opts.oversampling=1; end
    before=rng; rejected=false;
    try
        pod_basis({thin_data;excluded},[1;0],3,thin_opts);
    catch err
        rejected=contains(err.message,'target_rank exceeds the numerical rank');
    end
    result=test_true(result,['excluded rank ' method{1}],rejected&&isequal(before,rng),...
                     'unsupported retained rank fails without inventing modes or changing caller RNG');
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

% Prepare a sparse coordinate-selection model without expanding basis or state zeros
selected=pod_project(sparse_basis,{speye(2048)},sparse_state,1i*sparse_state);
result=test_true(result,'sparse prepared storage',issparse(selected.terms{1})&&...
                 issparse(selected.rho)&&issparse(selected.coil),...
                 'sparse coordinate selection preserves reduced term and state representations');
result=test_close(result,'sparse prepared term',selected.terms{1},eye(2),1e-12,0,...
                  'sparse preparation reproduces the independent identity projection');
result=test_close(result,'sparse prepared observable',pod_evolve(selected,1,0.03,0),...
                  -1i,0,0,'sparse complex detection retains the Hermitian initial overlap');
result=test_close(result,'sparse prepared evolution',pod_evolve(selected,1,0.03,3),...
                  -1i*exp(-1i*0.03*(0:3)),1e-12,0,...
                  'nonzero-duration sparse models recover the independent scalar phase trajectory');

% Reject nonfinite dense terms without extracting their double entries
for bad_value=[NaN Inf]
    bad_term=eye(3); bad_term(1,2)=bad_value; rejected=false;
    try
        pod_project(eye(3),{bad_term},ones(3,1),ones(3,1));
    catch err
        rejected=contains(err.message,'each term must be a finite');
    end
    result=test_true(result,['dense nonfinite ' num2str(bad_value)],rejected,...
                     'dense validation rejects NaN and Inf while preserving sparse support');
end

% Reject computed state overflow while admitting large finite projections
for storage={'dense','sparse'}
    subspace=[1;1]/sqrt(2); finite=[1;0];
    if strcmp(storage{1},'sparse'), subspace=sparse(subspace); finite=sparse(finite); end
    for phase=[1 1i]
        huge=phase*realmax*[1;1];
        if strcmp(storage{1},'sparse'), huge=sparse(huge); end
        for component={'rho','coil'}
            rho=finite; coil=finite; rejected=false;
            if strcmp(component{1},'rho'), rho=huge; else, coil=huge; end
            try
                pod_project(subspace,{eye(2)},rho,coil);
            catch err
                rejected=contains(err.message,'projected terms and states must be finite');
            end
            result=test_true(result,['projected overflow ' storage{1} ' ' component{1} ' ' num2str(phase)],...
                             rejected,'finite caller vectors whose reduced coordinates overflow fail during preparation');
        end
    end
    large=(realmax/2)*[1;1];
    if strcmp(storage{1},'sparse'), large=sparse(large); end
    prepared=pod_project(subspace,{eye(2)},large,1i*large);
    result=test_close(result,['finite large rho ' storage{1}],prepared.rho/realmax,...
                      1/sqrt(2),1e-15,0,'large finite reduced initial coordinates remain supported');
    result=test_close(result,['finite large coil ' storage{1}],prepared.coil/realmax,...
                      1i/sqrt(2),1e-15,0,'large finite complex detector coordinates remain supported');
end

% Compare streamed dense projection with independently selected coordinates
selection=[1;333;1000]; subspace=sparse(selection,(1:3).',1,1000,3);
dense_term=randn(1000)+1i*randn(1000); rho_dense=randn(1000,1); coil_dense=1i*randn(1000,1);
prepared=pod_project(subspace,{dense_term},rho_dense,coil_dense);
result=test_close(result,'streamed dense term',prepared.terms{1},dense_term(selection,selection),...
                  0,0,'coordinate selection gives an independent dense-operator projection reference');
result=test_close(result,'streamed dense state',prepared.rho,rho_dense(selection),...
                  0,0,'state projection uses the same bounded normalisation arithmetic');
result=test_close(result,'streamed dense detector',prepared.coil,coil_dense(selection),...
                  0,0,'Hermitian detector projection preserves the selected complex coordinates');
clear dense_term rho_dense coil_dense

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

% Check a zero-step request against the independently known initial overlap
result=test_close(result,'dense zero coefficients',pod_evolve(model,zeros(2,1),0.03,3),...
                  repmat(coil'*rho,1,4),1e-12,0,...
                  'exact zero coefficients repeat the initial complex overlap for positive time');
result=test_close(result,'zero-step initial observable',pod_evolve(model,coefficients,0.03,0),...
                  reference(1),1e-12,0,'zero steps return exactly one initial complex observable');

% Detect in normalised coordinates before restoring the physical observable
for storage={'dense','sparse'}
    subspace=eye(2); terms={zeros(2)}; rho=[2;-2];
    if strcmp(storage{1},'sparse')
        subspace=sparse(subspace); terms={sparse(2,2)}; rho=sparse(rho);
    end
    for phase=[1 1+1i]
        coil=phase*realmax*[1;1];
        if strcmp(storage{1},'sparse'), coil=sparse(coil); end
        prepared=pod_project(subspace,terms,rho,coil);
        for query={[0 0 0],[1 0 3],[0 0.03 3],[1 0.03 3]}
            evaluation=query{1};
            result=test_close(result,['scaled cancelling overlap ' storage{1} ' ' num2str(phase) ' ' num2str(evaluation)],...
                              pod_evolve(prepared,evaluation(1),evaluation(2),evaluation(3)),...
                              zeros(1,evaluation(3)+1),0,0,...
                              'bounded detection cancels large opposing products before physical rescaling');
        end
        prepared.rho=abs(rho); rejected=false;
        try
            pod_evolve(prepared,1,0.03,3);
        catch err
            rejected=contains(err.message,'detected observable must be finite');
        end
        result=test_true(result,['noncancelling overflow ' storage{1} ' ' num2str(phase)],...
                         rejected,'an unrepresentable physical observable remains an explicit error');
    end
end

% Zero-duration traces depend only on the initial overlap, not generator assembly
large_model=struct('terms',{{realmax*eye(2)}},'rho',[1;2],'coil',[1i;3]);
result=test_close(result,'zero-duration large generator',pod_evolve(large_model,2,0,3),...
                  repmat(6-1i,1,4),0,0,'a zero-duration trace is well defined even when affine assembly would overflow');

% Validate sparse reduced terms even when no generator is needed
zero_coeff_model=struct('terms',{{speye(100000)}},...
                        'rho',sparse(1,1,2+1i,100000,1),...
                        'coil',sparse(1,1,3+2i,100000,1));
result=test_close(result,'large sparse zero coefficients',...
                  pod_evolve(zero_coeff_model,0,0.03,3),repmat(8-1i,1,4),0,0,...
                  'zero coefficients avoid dense assembly while retaining complex conjugated detection');
sparse_prepared=struct('terms',{{speye(2048);sparse(2048,2048)}},...
                       'rho',ones(2048,1),'coil',1i*ones(2048,1));
result=test_close(result,'sparse zero-step observable',pod_evolve(sparse_prepared,[1;0],0.03,0),...
                  -2048i,0,0,'sparse reduced terms permit initial-only requests without expansion');
result=test_close(result,'sparse zero-duration trace',pod_evolve(sparse_prepared,[1;0],0,3),...
                  -2048i*ones(1,4),0,0,'sparse reduced terms and a zero term preserve repeated initial overlaps');
for bad_value=[NaN Inf]
    invalid_prepared=sparse_prepared;
    invalid_prepared.terms{1}=sparse(1,1,bad_value,2048,2048); rejected=false;
    try
        pod_evolve(invalid_prepared,[1;0],0,0);
    catch err
        rejected=contains(err.message,'model terms must be conformable finite');
    end
    result=test_true(result,['sparse model nonfinite ' num2str(bad_value)],rejected,...
                     'stored NaN and Inf are rejected before an initial-only return');
    rejected=false;
    try
        pod_evolve(invalid_prepared,[0;0],0.03,3);
    catch err
        rejected=contains(err.message,'model terms must be conformable finite');
    end
    result=test_true(result,['zero coefficients nonfinite ' num2str(bad_value)],rejected,...
                     'zero coefficients do not bypass validation of the supplied terms');
end

% Return exact constant signals without assembling unnecessary dense propagators
constant_model=struct('terms',{{speye(4);speye(4)}},...
                     'rho',[1+1i;2;0;0],'coil',[-1i;3;0;0]);
constant_trace=(5+1i)*ones(1,4);
result=test_close(result,'zero coefficients',pod_evolve(constant_model,[0;0],0.03,3),...
                  constant_trace,0,0,'all-zero coefficients give the initial observable at every time');
result=test_close(result,'sparse zero coefficients',pod_evolve(constant_model,sparse(2,1),0.03,3),...
                  constant_trace,0,0,'a sparse zero coefficient column has the same affine meaning');
zero_terms=constant_model; zero_terms.terms={sparse(4,4);sparse(4,4)};
result=test_close(result,'zero generator terms',pod_evolve(zero_terms,[2;-3],0.03,3),...
                  constant_trace,0,0,'zero sparse terms give a zero generator despite nonzero coefficients');
result=test_close(result,'cancelled generator',pod_evolve(constant_model,[1;-1],0.03,3),...
                  constant_trace,0,0,'exact affine cancellation gives a constant signal');
zero_state=constant_model; zero_state.rho=zeros(4,1);
result=test_close(result,'zero initial state',pod_evolve(zero_state,[1;1],0.03,3),...
                  zeros(1,4),0,0,'linear propagation preserves an exactly zero initial state');
zero_detector=constant_model; zero_detector.coil=sparse(4,1);
result=test_close(result,'zero detector',pod_evolve(zero_detector,[1;1],0.03,3),...
                  zeros(1,4),0,0,'an exactly zero detector gives no signal for any finite generator');
mixing_model=struct('terms',{{[0 1;1 0]}},'rho',[1;0],'coil',[0;1]);
result=test_close(result,'initial zero is not constant',pod_evolve(mixing_model,1,0.1,3),...
                  -1i*sin(0.1*(0:3)),1e-12,0,...
                  'zero initial overlap alone does not imply a zero trajectory under state mixing');

% Check exact sparse right and left nullspaces without requiring a zero generator
right_null=struct('terms',{{sparse([0 4 0;0 0 3;0 0 2])}},...
                  'rho',[1;0;0],'coil',[1+1i;2;0]);
left_null=struct('terms',{{sparse([0 0 0;4 0 3;0 0 2])}},...
                 'rho',[1+1i;1;1],'coil',[1;0;0]);
result=test_close(result,'exact stationary state',pod_evolve(right_null,1,0.03,3),...
                  (1-1i)*ones(1,4),0,0,'a nonzero generator leaves its exact right-nullspace state unchanged');
result=test_close(result,'exact stationary detector',pod_evolve(left_null,1,0.03,3),...
                  (1+1i)*ones(1,4),0,0,'a left-nullspace detector is constant despite nonstationary state evolution');
small_phase=struct('terms',{{sparse(1e-100)}},'rho',1,'coil',1);
small_fid=pod_evolve(small_phase,1,0.03,3);
result=test_true(result,'nonzero small generator',imag(small_fid(end))~=0,...
                 'stationary returns use structural support rather than a magnitude tolerance');

% Reject vanished normalised products without confusing exact cancellation with loss
for storage={'dense','sparse'}
    for phase=[1 1i]
        tiny_basis=[1;phase*2^-1074]; huge_state=[0;2^1023]; huge_term=[0 2^1023;0 0];
        if strcmp(storage{1},'sparse')
            tiny_basis=sparse(tiny_basis); huge_state=sparse(huge_state); huge_term=sparse(huge_term);
        end
        for query={@()pod_project(tiny_basis,{zeros(2)},huge_state,[1;0]),...
                   @()pod_project(tiny_basis,{zeros(2)},[1;0],huge_state),...
                   @()pod_project(tiny_basis,{huge_term},[1;0],[1;0])}
            rejected=false;
            try
                query{1}();
            catch err
                rejected=contains(err.message,'product dynamic range');
            end
            result=test_true(result,['projection product loss ' storage{1} ' ' num2str(phase)],...
                             rejected,'tiny basis factors lost inside either matrix product are rejected');
        end
        prepared=struct('terms',{{zeros(2)}},'rho',[2^1023;phase*2^-50],'coil',[0;2^1023]);
        if strcmp(storage{1},'sparse')
            prepared.rho=sparse(prepared.rho); prepared.coil=sparse(prepared.coil);
        end
        for query={[0 0 0],[1 0 2],[0 0.03 2],[1 0.03 2]}
            evaluation=query{1};
            result=test_close(result,['detection product preservation ' storage{1} ' ' num2str(phase)],...
                              pod_evolve(prepared,evaluation(1),evaluation(2),evaluation(3))/2^973,...
                              phase*ones(1,evaluation(3)+1),0,0,...
                              'exponent-aware detection retains the finite overlap in every constant query');
        end
    end
    prepared=struct('terms',{{0.5+1i*2^-1074}},'rho',1,'coil',1);
    if strcmp(storage{1},'sparse'), prepared.terms={sparse(prepared.terms{1})}; end
    rejected=false;
    try
        pod_evolve(prepared,1i,1,1);
    catch err
        rejected=contains(err.message,'product dynamic range');
    end
    result=test_true(result,['affine component product loss ' storage{1}],rejected,...
                     'a surviving complex component cannot mask a vanished constituent product');
    prepared=struct('terms',{{0.5+0.5i}},'rho',1,'coil',1);
    result=test_close(result,['complex component cancellation ' storage{1}],...
                      pod_evolve(prepared,1+1i,1,2),exp(0:2),1e-12,0,...
                      'representable real products may cancel exactly during complex multiplication');
    prepared=struct('terms',{{zeros(2)}},'rho',[1;1],'coil',[1;-1]);
    result=test_close(result,['detection exact cancellation ' storage{1}],...
                      pod_evolve(prepared,0,0,0),0,0,0,...
                      'a zero overlap from finite cancelling products remains supported');
    exact_basis=[1;1]/sqrt(2);
    if strcmp(storage{1},'sparse'), exact_basis=sparse(exact_basis); end
    prepared=pod_project(exact_basis,{zeros(2)},2^1000*[1;-1],[1;0]);
    result=test_close(result,['projection exact cancellation ' storage{1}],prepared.rho,0,0,0,...
                      'a zero projected coordinate from finite cancelling products remains supported');
end

% Preserve representable time-scaled phases across scalar factor ranges
for storage={'dense','sparse'}
    terms={1e-200;2^1000;-1i*2^-1000;-1i*2^800;2^-700;1i*700;1i*log(2)};
    if strcmp(storage{1},'sparse'), terms=cellfun(@sparse,terms,'UniformOutput',false); end
    prepared=pod_project(1,terms(1),1,1);
    observed=pod_evolve(prepared,1e-200,1e200,1);
    result=test_close(result,['affine tiny phase ' storage{1}],imag(observed(2))/1e-200,...
                      -1,1e-12,0,'the time step rescues a phase lost by coefficient-term multiplication');
    prepared.terms=terms(2); observed=pod_evolve(prepared,2^-700,2^-700,1);
    result=test_close(result,['scalar scale underflow ' storage{1}],imag(observed(2))/2^-400,...
                      -1,1e-12,0,'the binary frame avoids an underflowing coefficient-time product');
    prepared.terms=terms(3);
    result=test_close(result,['scalar scale overflow ' storage{1}],pod_evolve(prepared,2^700,2^700,1),...
                      [1 0],0,0,'the binary frame retains a finite strongly damped exponent');
    prepared.terms=terms(4);
    result=test_close(result,['generator scale overflow ' storage{1}],pod_evolve(prepared,2^400,2^-1000,1),...
                      [1 0],0,0,'the common binary frame avoids an overflowing unscaled generator');
    prepared.terms=terms(5); rejected=false;
    try
        pod_evolve(prepared,2^-700,2^-700,1);
    catch err
        rejected=contains(err.message,'represent');
    end
    result=test_true(result,['unrepresentable scaled term ' storage{1}],rejected,...
                     'lost nonzero entries are rejected rather than certified as structural zeros');

    % Preserve detector-scaled amplification and reject unrepresentable signals
    prepared.terms=terms(6); prepared.coil=1e-300;
    observed=pod_evolve(prepared,1,1,2);
    result=test_close(result,['scaled amplified state ' storage{1}],...
                      observed(3)/exp(1400+log(1e-300)),1,1e-12,0,...
                      'binary state amplitude retains a finite observable beyond physical state range');
    observed=pod_evolve(prepared,1,1,1);
    result=test_close(result,['finite amplified state ' storage{1}],observed(2)/exp(700+log(1e-300)),...
                      1,1e-12,0,'the preceding representable state and observable remain supported');
    prepared.terms=terms(7);
    for phase=[1 1+1i]
        prepared.coil=phase*(realmax/2); rejected=false;
        try
            pod_evolve(prepared,1,1,2);
        catch err
            rejected=contains(err.message,'detected observable must be finite');
        end
        result=test_true(result,['later detector overflow ' storage{1} ' ' num2str(phase)],rejected,...
                         'finite propagated states whose detection overflows are rejected');
        prepared.coil=phase*(realmax/4);
        result=test_close(result,['finite amplified detector ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,1,1,2)/realmax,conj(phase)*[1/4 1/2 1],1e-12,0,...
                          'finite real and complex amplified observables remain supported');
    end
end

% Retain finite factor orderings and cancellation at the exponent-range boundary
for storage={'dense','sparse'}
    terms={2^1023;realmax;diag([0 realmax])};
    if strcmp(storage{1},'sparse'), terms=cellfun(@sparse,terms,'UniformOutput',false); end
    prepared=pod_project(1,terms(1),1,1);
    result=test_close(result,['subnormal scale ordering ' storage{1}],...
                      pod_evolve(prepared,2^51,2^-1074,1),[1 exp(-1i)],1e-12,0,...
                      'the common binary frame retains a finite exponent despite coefficient-term overflow');
    observed=pod_evolve(prepared,2.1,2^-1074,1);
    result=test_close(result,['subnormal scale precision ' storage{1}],imag(observed(2))/(2.1*2^-51),...
                      -1,1e-12,0,'subnormal scalar products do not quantise a representable phase');
    result=test_close(result,['subnormal coefficient ordering ' storage{1}],...
                      pod_evolve(prepared,2^-1074,2^40,1),[1 exp(-1i*2^-11)],1e-12,0,...
                      'separate exponents retain the phase without overflowing time-term multiplication');
    prepared.terms=terms([2;2]);
    result=test_close(result,['overflowing scaled cancellation ' storage{1}],...
                      pod_evolve(prepared,[1;-1],2,3),ones(1,4),0,0,...
                      'finite affine cancellation precedes an otherwise overflowing time scale');
    prepared=pod_project(eye(2),terms(3),[1;0],[1;0]);
    result=test_close(result,['unreachable state overflow ' storage{1}],...
                      pod_evolve(prepared,2,1,3),ones(1,4),0,0,...
                      'common structural stationarity precedes irrelevant-entry scaling');
    prepared.rho=[0;1];
    result=test_close(result,['unreachable detector overflow ' storage{1}],...
                      pod_evolve(prepared,2,1,3),zeros(1,4),0,0,...
                      'a common stationary detector needs no unreachable propagator entries');
end

% Require the same stationary side across all terms, not a different side per term
mixed_sides=pod_project(eye(3),{sparse(3,2,1,3,3);sparse(2,1,1,3,3)},[1;0;0],[0;0;1]);
result=test_close(result,'mixed stationary sides',pod_evolve(mixed_sides,[1;1],0.03,3),...
                  -(0.03*(0:3)).^2/2,1e-12,0,...
                  'two individually invisible terms compose a detector-visible nilpotent chain');

% Retain a detector-visible tiny phase of a sparse mixing generator
scaled_mixing=pod_project(eye(2),{1e-200*sparse([0 1;1 0])},[1;0],[0;1]);
observed=pod_evolve(scaled_mixing,1e-200,1e200,1);
result=test_close(result,'affine tiny mixing phase',imag(observed(2))/1e-200,-1,1e-12,0,...
                  'a zero initial overlap does not hide a representable off-diagonal phase');

% Check binary-frame scaling and common structural stationarity
for storage={'dense','sparse'}
    for phase=[1 1i]
        huge=phase*2^1023; stationary=phase*[0 realmax;0 realmax];
        if strcmp(storage{1},'sparse'), huge=sparse(huge); stationary=sparse(stationary); end
        prepared=struct('terms',{{huge}},'rho',1,'coil',1);
        result=test_close(result,['subnormal fallback ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,2^51,2^-1074,2),exp(-1i*phase*(0:2)),1e-12,0,...
                          'a common binary frame avoids overflowing term-first multiplication');
        prepared.terms={phase*realmax;-phase*realmax;phase*0.5};
        if strcmp(storage{1},'sparse'), prepared.terms=cellfun(@sparse,prepared.terms,'UniformOutput',false); end
        result=test_close(result,['affine scaling fallback ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,[1;1;1],2,2),exp(-1i*phase*(0:2)),1e-12,0,...
                          'unscaled cancellation retains a finite nonzero exponent when separate scaling overflows');
        prepared=struct('terms',{{stationary;stationary}},'rho',[1;0],'coil',[1;1]);
        result=test_close(result,['unreachable state entries ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,[2;3],1,2),ones(1,3),0,0,...
                          'common zero state-support columns avoid irrelevant range overflow');
        prepared.terms={stationary.';stationary.'}; prepared.rho=[1;1]; prepared.coil=[1;0];
        result=test_close(result,['unreachable detector entries ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,[2;3],1,2),ones(1,3),0,0,...
                          'common zero detector-support rows avoid irrelevant range overflow');
    end
    prepared=struct('terms',{{[0 1;0 0];[0 0;1 0]}},'rho',[1;0],'coil',[1;0]);
    if strcmp(storage{1},'sparse'), prepared.terms=cellfun(@sparse,prepared.terms,'UniformOutput',false); end
    result=test_close(result,['mixed stationarity sides ' storage{1}],...
                      pod_evolve(prepared,[1;1],0.2,2),cos(0.2*(0:2)),1e-12,0,...
                      'state stationarity in one term and detector stationarity in another do not combine');
    prepared.terms={realmax;-realmax;2^-1000}; prepared.rho=1; prepared.coil=1;
    if strcmp(storage{1},'sparse'), prepared.terms=cellfun(@sparse,prepared.terms,'UniformOutput',false); end
    rejected=false;
    try
        pod_evolve(prepared,[1;1;1],2^-1000,1);
    catch err
        rejected=contains(err.message,'represent');
    end
    result=test_true(result,['fallback underflow rejected ' storage{1}],rejected,...
                     'lost nonzero term entries remain errors rather than stationary signals');
end

% Preserve finite phases when a state action underflows rather than vanishes
scaled_state=struct('terms',{{sparse(1e-200)}},'rho',1e-200,'coil',1e200);
result=test_close(result,'underflowing state action',pod_evolve(scaled_state,1,1e200,3),...
                  exp(-1i*(0:3)),1e-12,0,...
                  'an underflowing action is not a structurally stationary state');
scaled_coil=struct('terms',{{sparse(1e-200)}},'rho',1e200,'coil',1e-200);
result=test_close(result,'underflowing detector action',pod_evolve(scaled_coil,1,1e200,3),...
                  exp(-1i*(0:3)),1e-12,0,...
                  'an underflowing action is not a structurally stationary detector');
cancel_state=struct('terms',{{sparse([1 -1;2 -2])}},'rho',[1;1],'coil',[1;0]);
result=test_close(result,'cancellation stationary state',pod_evolve(cancel_state,1,0.03,3),...
                  ones(1,4),1e-12,0,...
                  'a nullspace from cancellation remains supported through dense propagation');

% Compare amplified non-normal composition with a scaled analytic reference
amplified=pod_project(eye(2),{1i*[700 1;0 699]},[0;1],[1e-300;0]);
observed=pod_evolve(amplified,1,1,2); times=0:2;
reference=exp(700*times+log(1e-300)).*(-expm1(-times));
result=test_close(result,'amplified non-normal composition',observed(2:3)./reference(2:3),...
                  ones(1,2),1e-12,0,...
                  'binary amplitude propagation preserves a detector-scaled triangular semigroup');

% Guard scaled matrix products and retain exponent-aware observable products
for storage={'dense','sparse'}
    subspace=[1;2^-1074]; initial=[0;2^1023]; detector=[1;0]; terms={zeros(2)};
    if strcmp(storage{1},'sparse')
        subspace=sparse(subspace); initial=sparse(initial); detector=sparse(detector); terms={sparse(2,2)};
    end
    rejected=false;
    try
        pod_project(subspace,terms,initial,detector);
    catch err
        rejected=contains(err.message,'product dynamic range');
    end
    result=test_true(result,['projection product range ' storage{1}],rejected,...
                     'a tiny basis component must not be silently lost in a scaled product');
    prepared=struct('terms',{{0.5+1i*2^-1074}},'rho',1,'coil',1);
    if strcmp(storage{1},'sparse'), prepared.terms={sparse(prepared.terms{1})}; end
    rejected=false;
    try
        pod_evolve(prepared,1i,1,1);
    catch err
        rejected=contains(err.message,'product dynamic range');
    end
    result=test_true(result,['affine component product range ' storage{1}],rejected,...
                     'nonzero complex components must not disappear before a range check');
    for phase=[1 1i 1+1i]
        prepared=struct('terms',{terms},'rho',[2^1023;2^-50],'coil',phase*[0;2^1023]);
        if strcmp(storage{1},'sparse')
            prepared.rho=sparse(prepared.rho); prepared.coil=sparse(prepared.coil);
        end
        result=test_close(result,['exponent-aware overlap ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,0,0,0)/2^973,conj(phase),0,0,...
                          'component exponents preserve a finite dot product across vector scales');
    end
end

% Remove inaccessible growth while retaining directed paths to the observable
for storage={'dense','sparse'}
    terms={1i*diag([2000 -1])}; subspace=eye(2);
    if strcmp(storage{1},'sparse'), terms={sparse(terms{1})}; subspace=sparse(subspace); end
    prepared=pod_project(subspace,terms,[0;1],[0;1]);
    result=test_close(result,['disconnected growth ' storage{1}],pod_evolve(prepared,1,1,3),...
                      exp(-(0:3)),1e-14,0,'an unreachable unstable block cannot spoil a finite observable');
    A=diag([-1 -2 -3 2000]); A(1,2)=1; A(2,3)=1;
    if strcmp(storage{1},'sparse'), A=sparse(A); end
    prepared=struct('terms',{{1i*A}},'rho',[0;0;1;0],'coil',[1;0;0;0]);
    times=0:3; reference=exp(-times).*(-expm1(-times)).^2/2;
    result=test_close(result,['directed observable paths ' storage{1}],pod_evolve(prepared,1,1,3),...
                      reference,1e-13,0,'all intermediate coordinates on a source-to-detector path remain');
end

% Compare dense invariant actions with exact modes and preserve weak growing controls
for storage={'dense','sparse'}
    A=[999.5 1000.5;1000.5 999.5]; subspace=eye(2); before=rng;
    if strcmp(storage{1},'sparse'), A=sparse(A); subspace=sparse(subspace); end
    for phase=[1 1i]
        prepared=pod_project(subspace,{1i*A},phase*[1;-1],[1;-1]);
        for query={[1 3],[2^-10 1024]}
            evaluation=query{1}; times=(0:evaluation(2))*evaluation(1);
            result=test_close(result,['dense invariant action ' storage{1} ' ' num2str(phase)],...
                              pod_evolve(prepared,1,evaluation(1),evaluation(2)),...
                              2*phase*exp(-times),1e-11,0,...
                              'a dense unobserved growth mode is never exponentiated as a full propagator');
        end
        prepared=pod_project(subspace,{1i*(A+1i*[1.5 -0.5;-0.5 1.5])},[1;0],phase*[1;-1]);
        result=test_close(result,['complex adjoint action ' storage{1} ' ' num2str(phase)],...
                          pod_evolve(prepared,1,1,3),conj(phase)*exp((-1+2i)*(0:3)),1e-11,0,...
                          'adjoint detection retains the propagation sign and Hermitian conjugation');
    end
    result=test_true(result,['native action RNG ' storage{1}],isequal(before,rng),...
                     'native norm estimation must not consume the caller random stream');
    prepared=pod_project(subspace,{1i*A},[1;1],[1;1]); rejected=false; before=rng;
    try
        pod_evolve(prepared,1,1,1);
    catch err
        rejected=contains(err.message,'detected observable must be finite');
    end
    result=test_true(result,['observed growing mode ' storage{1}],rejected,...
                     'a genuinely unrepresentable growing observable remains an error');
    result=test_true(result,['failed action RNG ' storage{1}],isequal(before,rng),...
                     'native norm-estimation errors also preserve the caller random stream');
    prepared.rho=[1+eps;-1+eps]; prepared.coil=prepared.rho; rejected=false;
    try
        pod_evolve(prepared,1,1,1);
    catch
        rejected=true;
    end
    result=test_true(result,['weak growing mode control ' storage{1}],rejected,...
                     'an arbitrarily small representable growing contribution is not dropped as invariant');
end

% Compare complete-support normal dynamics with an independent rank-one formula
for storage={'dense','sparse'}
    L=-ones(3)-1i*eye(3); subspace=eye(3); rho=[1;0;0]; times=(0:8)*0.03;
    if strcmp(storage{1},'sparse'), L=sparse(L); subspace=sparse(subspace); end
    for phase=[1 1i]
        coil=phase*[1;2;-1];
        reference=conj(phase)*exp(-times).*(1+2*(exp(3*1i*times)-1)/3);
        for n=1:2
            if n==1, terms={L}; coefficients=1;
            else, terms={eye(3);L-eye(3)}; coefficients=[1;1]; end
            prepared=pod_project(subspace,terms,rho,coil);
            result=test_close(result,['complete support action ' storage{1} ' ' num2str(phase) ' ' num2str(n)],...
                              pod_evolve(prepared,coefficients,0.03,8),reference,1e-12,0,...
                              'complete term support and cached component bounds preserve the normal semigroup');
        end
    end
end

% Compare dense-basis diagonal projection with independent weighted parity sums
for storage={'dense','sparse'}
    nranks=128; diagonal=(1:nranks)'/nranks; term=diag(diagonal);
    if strcmp(storage{1},'sparse'), term=sparse(term); end
    for phase=[1 1i]
        subspace=[ones(nranks,1),phase*((-1).^(1:nranks)')]/sqrt(nranks);
        prepared=pod_project(subspace,{term},eye(nranks,1),1i*eye(nranks,1));
        average=(nranks+1)/(2*nranks); coordinates=[1;-conj(phase)]/sqrt(nranks);
        reference=[average phase/(2*nranks);conj(phase)/(2*nranks) average];
        result=test_close(result,['fixed basis range ' storage{1} ' ' num2str(phase)],...
                          [prepared.terms{1},prepared.rho,prepared.coil],...
                          [reference,coordinates,1i*coordinates],1e-12,0,...
                          'dense fixed coordinates preserve independent diagonal and parity projections');
    end
end

% Compare incomplete dense support and disconnected controls with full exponentials
for storage={'dense','sparse'}
    for removed=[1 4]
        L=reshape((1:25)/25,5,5)+1i*eye(5); L(2:removed+1,1)=0; subspace=eye(5);
        if strcmp(storage{1},'sparse'), L=sparse(L); subspace=sparse(subspace); end
        for phase=[1 1i]
            rho=eye(5,1); coil=phase*[0;1;0;0;0]; reference=complex(zeros(1,5));
            for n=0:4, reference(n+1)=coil'*expm(-1i*full(L)*(0.03*n))*rho; end
            prepared=pod_project(subspace,{L},rho,coil);
            result=test_close(result,['implicit incomplete support ' storage{1} ' ' num2str(removed) ' ' num2str(phase)],...
                              pod_evolve(prepared,1,0.03,4),reference,1e-12,0,...
                              'implicit directed traversal retains paths without assuming dense connectivity');
        end
    end
end

% Preserve real projection storage without discarding genuinely complex coordinates
for basis_kind={'dense','sparse'}
    for term_kind={'dense','sparse'}
        subspace=eye(3,2); term=diag(1:3); rho=[1;2;3]; coil=[3;2;1];
        if strcmp(basis_kind{1},'sparse'), subspace=sparse(subspace); end
        if strcmp(term_kind{1},'sparse'), term=sparse(term); end
        label=[basis_kind{1} ' ' term_kind{1}];
        prepared=pod_project(subspace,{term},rho,coil);
        result=test_true(result,['real projected storage ' label],...
                         isreal(prepared.terms{1})&&isreal(prepared.rho)&&isreal(prepared.coil),...
                         'all-real Galerkin data retain real prepared arrays');
        result=test_close(result,['real projected values ' label],...
                          [prepared.terms{1},prepared.rho,prepared.coil],...
                          [1 0 1 3;0 2 2 2],0,0,...
                          'coordinate selection agrees with independent real entries');
        term(1,2)=1i; prepared=pod_project(subspace,{term},1i*rho,coil);
        result=test_close(result,['complex storage control ' label],...
                          [prepared.terms{1},prepared.rho,prepared.coil],...
                          [1 1i 1i 3;0 2 2i 2],0,0,...
                          'complex term and state contributions survive real initial allocation');
    end
end

% Preserve scalar phase independently of rounded Rayleigh-rate selection
for frequency=[0.3 0.7 -0.3 -0.7]
    for swapped=[false true]
        rho=0.3+0.4i; coil=-0.4+0.7i;
        if swapped, [rho,coil]=deal(coil,rho); end
        prepared=pod_project(1,{frequency},rho,coil);
        reference=(coil'*rho)*exp(-1i*frequency*(0:256)*0.01);
        result=test_close(result,['scalar phase convention ' num2str(frequency) ' ' num2str(swapped)],...
                          pod_evolve(prepared,1,0.01,256),reference,1e-12,0,...
                          'the closed-form scalar action preserves phase and Hermitian detection');
    end
end

% Compare real relaxation actions with the independent rank-one exponential
for storage={'dense','sparse'}
    ranks=4; generator=-eye(ranks)+ones(ranks)/(2*ranks); subspace=eye(ranks);
    if strcmp(storage{1},'sparse'), generator=sparse(generator); subspace=sparse(subspace); end
    for coefficient=[1 1+1i]
        for phase=[1 1i]
            prepared=pod_project(subspace,{1i*generator},ones(ranks,1),phase*ones(ranks,1));
            reference=conj(phase)*ranks*exp(-coefficient*(0:8)*0.1/2);
            result=test_close(result,['real relaxation action ' storage{1} ' ' num2str(coefficient) ' ' num2str(phase)],...
                              pod_evolve(prepared,coefficient,0.1,8),reference,1e-12,0,...
                              'real and genuinely complex actions preserve the rank-one relaxation semigroup');
        end
    end
end

% Compare separated complex diagonals and adjoints with a triangular closed form
for storage={'dense','sparse'}
    generator=[-2+0.3i 3+2i;0 -1+0.5i]; subspace=eye(2); times=(0:8)*0.03;
    first=exp(generator(1,1)*times); last=exp(generator(2,2)*times);
    transfer=generator(1,2)*(first-last)/(generator(1,1)-generator(2,2));
    if strcmp(storage{1},'sparse'), generator=sparse(generator); subspace=sparse(subspace); end
    for swapped=[false true]
        rho=[0.5;1]; coil=[1;0.3i];
        if swapped, [rho,coil]=deal(coil,rho); end
        reference=conj(coil(1))*(first*rho(1)+transfer*rho(2))+conj(coil(2))*last*rho(2);
        prepared=pod_project(subspace,{1i*generator},rho,coil);
        result=test_close(result,['separated diagonal action ' storage{1} ' ' num2str(swapped)],...
                          pod_evolve(prepared,1,0.03,8),reference,1e-12,0,...
                          'complex diagonal and off-diagonal contributions retain Hermitian detection');
    end
end

% Compare connected tridiagonal actions with an independent alternating eigenmode
for storage={'dense','sparse'}
    ranks=128;
    if strcmp(storage{1},'sparse'), ranks=8192; end
    values=ones(ranks,1); diagonal=-2*values; diagonal([1 end])=-3;
    generator=spdiags([values diagonal values],-1:1,ranks,ranks);
    if strcmp(storage{1},'dense'), generator=full(generator); end
    for coefficient=[1 1+0.2i]
        for swapped=[false true]
            rho=(-1).^(1:ranks).'; coil=(1+0.3i)*rho;
            if swapped, [rho,coil]=deal(coil,rho); end
            prepared=struct('terms',{{1i*generator}},'rho',rho,'coil',coil);
            reference=(coil'*rho)*exp(-4*coefficient*(0:4)*0.01);
            result=test_close(result,['tridiagonal action ' storage{1} ' ' num2str(coefficient) ' ' num2str(swapped)],...
                              pod_evolve(prepared,coefficient,0.01,4),reference,0,1e-12,...
                              'connected sparse and dense actions preserve the exact relaxation eigenmode');
        end
    end
end

% Remove disconnected affine cancellations before restoring irrelevant growth
for storage={'dense','sparse'}
    for residual=[0 1]
        for phase=[1 1i]
            first=[1 1;1 1i*2^1023]; second=[0 -1;-1 residual*1i*2^1023];
            if strcmp(storage{1},'sparse'), first=sparse(first); second=sparse(second); end
            prepared=struct('terms',{{first;second}},'rho',[1;0],'coil',[phase;0]);
            reference=conj(phase)*exp(-1i*(0:3));
            result=test_close(result,['affine disconnected ' storage{1} ' ' num2str(residual) ' ' num2str(phase)],...
                              pod_evolve(prepared,[1;1],1,3),reference,1e-12,0,...
                              'cancelled connecting paths exclude finite or unrepresentable irrelevant growth');
        end
    end
    second=[0 -0.5;-0.5 0];
    if strcmp(storage{1},'sparse'), second=sparse(second); end
    prepared.terms={first;second}; rejected=false;
    try
        pod_evolve(prepared,[1;1],1,3);
    catch
        rejected=true;
    end
    result=test_true(result,['affine connected range ' storage{1}],rejected,...
                     'a genuine path to extreme growth remains outside the supported working range');
end

% Preserve weak affine couplings between cancelling order-one contributions
for storage={'dense','sparse'}
    terms={[1 0;2^-54 0];[0 0;1 0];[-1 0;-1 0]}; orders=perms(1:3);
    if strcmp(storage{1},'sparse'), terms=cellfun(@sparse,terms,'UniformOutput',false); end
    for k=1:size(orders,1)
        for phase=[1 1i]
            prepared=struct('terms',{terms(orders(k,:))},'rho',[phase;0],'coil',[0;1+0.3i]);
            reference=(1-0.3i)*phase*(-1i)*(0:3);
            result=test_close(result,['affine residual ' storage{1} ' ' num2str(k) ' ' num2str(phase)],...
                              pod_evolve(prepared,ones(3,1),2^54,3),reference,1e-12,0,...
                              'all term permutations preserve the independently known nilpotent coupling');
        end
    end
end

% Reject uncertified zeros after a second-level compensated residual is lost
for storage={'dense','sparse'}
    terms={[1 0;1 0];[1 0;2^-54 0];[-1 0;2^-108 1];[-1 0;-1 0];[0 0;-2^-54 -1]};
    if strcmp(storage{1},'sparse'), terms=cellfun(@sparse,terms,'UniformOutput',false); end
    prepared=struct('terms',{terms},'rho',[1;0],'coil',[0;1]); rejected=false;
    try
        pod_evolve(prepared,ones(5,1),2^108,3);
    catch caught
        rejected=contains(caught.message,'affine cancellation cannot certify');
    end
    result=test_true(result,['uncertified affine zero ' storage{1}],rejected,...
                     'a rounded zero with a lost correction is an explicit unsupported sum');
    prepared.terms=terms([1 4 2 5 3]);
    result=test_close(result,['certified affine control ' storage{1}],...
                      pod_evolve(prepared,ones(5,1),2^108,3),-1i*(0:3),1e-12,0,...
                      'the independently known nilpotent signal survives an exactly accumulated control');
end

% Reject physical restoration loss without discarding representable subnormals
for storage={'dense','sparse'}
    for phase=[1 1i]
        subspace=[1;2^-1073]; initial=[1;0];
        term=[0 0;phase*realmin 0]; zero=zeros(2); tiny_state=[0;phase*realmin];
        if strcmp(storage{1},'sparse')
            subspace=sparse(subspace); term=sparse(term); zero=sparse(zero); tiny_state=sparse(tiny_state);
        end
        requests={@()pod_project(subspace,{term},initial,initial),...
                  @()pod_project(subspace,{zero},tiny_state,initial),...
                  @()pod_project(subspace,{zero},initial,tiny_state)};
        for n=1:numel(requests)
            rejected=false;
            try
                requests{n}();
            catch caught
                rejected=contains(caught.message,'physical double representation');
            end
            result=test_true(result,['restoration loss ' storage{1} ' ' num2str(phase) ' ' num2str(n)],rejected,...
                             'terms, initial coordinates, and detector coordinates cannot silently round to zero');
        end
        term=[phase*realmin 0;1i*phase*realmin 0];
        if strcmp(storage{1},'sparse'), term=sparse(term); end
        rejected=false;
        try
            pod_project(subspace,{term},initial,initial);
        catch caught
            rejected=contains(caught.message,'physical double representation');
        end
        result=test_true(result,['restoration component loss ' storage{1} ' ' num2str(phase)],rejected,...
                         'losing one complex component is rejected even when the coordinate remains nonzero');
        subspace=[1;2^-52]; term=[0 0;phase*realmin 0];
        if strcmp(storage{1},'sparse'), subspace=sparse(subspace); term=sparse(term); end
        prepared=pod_project(subspace,{term},initial,initial);
        result=test_true(result,['representable projected subnormal ' storage{1} ' ' num2str(phase)],...
                         isequal(full(prepared.terms{1}),phase*2^-1074),...
                         'the independent Galerkin coordinate at the smallest subnormal is preserved');
        if phase==1
            observed=pod_evolve(prepared,1,realmax,1); expected=imag(exp(-1i*(2^-1074*realmax)));
            result=test_close(result,['subnormal projected phase ' storage{1}],...
                              imag(observed(2))/expected,1,1e-12,0,...
                              'a long acquisition preserves the independently representable phase, not just its real background');
        end
    end
end

% Compare componentwise diagonal phases with independent scalar exponentials
for storage={'dense','sparse'}
    for frequency=[1 1e16]
        for swapped=[false true]
            for coefficient=[1 1+1i*2^-54]
                term=diag([frequency -frequency 0]); initial=[1;2;1i]; detector=[2;1i;1];
                if swapped, [initial,detector]=deal(detector,initial); end
                if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(3),{term},initial,detector);
                expected=complex(zeros(1,5));
                for k=0:4, expected(k+1)=detector'*(exp(-1i*coefficient*diag(term)*k).*initial); end
                saved_rng=rng; observed=pod_evolve(prepared,coefficient,1,4);
                result=test_close(result,['diagonal phases ' storage{1} ' ' num2str(frequency) ' ' num2str(swapped) ' ' num2str(coefficient)],...
                                  observed,expected,1e-12,0,...
                                  'diagonal modes use independent scalar exponentials without a generic action-norm limit');
                result=test_true(result,['diagonal RNG ' storage{1} ' ' num2str(frequency) ' ' num2str(swapped) ' ' num2str(coefficient)],...
                                 isequal(rng,saved_rng),'closed-form diagonal actions preserve the caller stream');
            end
        end
    end
end

% Scale each diagonal exponential jointly with its state and detector amplitudes
for storage={'dense','sparse'}
    term=1i*diag([1000 0]);
    if strcmp(storage{1},'sparse'), term=sparse(term); end
    prepared=pod_project(eye(2),{term},[1e-300;1],[1;1]);
    observed=pod_evolve(prepared,1,1,1); expected=1+exp(1000+log(1e-300));
    result=test_close(result,['diagonal finite product ' storage{1}],observed(2)/expected,1,1e-12,0,...
                      'a finite state-exponential product is not rejected because the exponential factor overflows');
    for rate=[-1000 -400 400 1000]
        amplitude=1; sensitivity=1;
        if rate>0, amplitude=1e-300; else, sensitivity=1e300; end
        nsteps=3; if rate==1000, nsteps=1; end
        for frequency=[0 1e16]
            for swapped=[false true]
                term=diag([frequency+1i*rate 0 0]);
                initial=[amplitude*(1+1i);1;1]; detector=[sensitivity;1;-1]; phase=1+1i;
                if swapped, [initial,detector]=deal(detector,initial); phase=conj(phase); end
                if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(3),{term},initial,detector);
                expected=phase*exp(log(amplitude)+log(sensitivity)+rate*(0:nsteps)).*...
                         exp(-1i*frequency*(0:nsteps));
                observed=pod_evolve(prepared,1,1,nsteps); nonzero=(expected~=0);
                observed(nonzero)=observed(nonzero)./expected(nonzero);
                result=test_close(result,['diagonal scaled modes ' storage{1} ' ' num2str(rate) ' ' num2str(frequency) ' ' num2str(swapped)],...
                                  observed,double(nonzero),1e-12,0,...
                                  'independent mode amplitudes retain tiny and large observables after exact background cancellation');
            end
        end
    end
end

% Certify physical damping tails without dropping representable collective tails
for storage={'dense','sparse'}
    for rate=[-750 -1000]
        term=1i*diag([rate 0]);
        if strcmp(storage{1},'sparse'), term=sparse(term); end
        for phase=[1 1+1i]
            prepared=pod_project(eye(2),{term},[1;1],[conj(phase);1]);
            result=test_close(result,['physical diagonal tail ' storage{1} ' ' num2str(rate) ' ' num2str(phase)],...
                              pod_evolve(prepared,1,1,3),1+phase*exp(rate*(0:3)),1e-12,0,...
                              'a certified physically unrepresentable damping tail does not hide the surviving signal');
        end
    end
    for count=[2 256]
        rate=-700; if count==256, rate=-750; end
        term=1i*diag([0;rate*ones(count,1)]); detector=[realmax/4;ones(count,1)];
        if strcmp(storage{1},'sparse'), term=sparse(term); end
        prepared=pod_project(eye(count+1),{term},ones(count+1,1),detector);
        rejected=false;
        try
            pod_evolve(prepared,1,1,1);
        catch caught
            rejected=contains(caught.message,'observable product dynamic range');
        end
        result=test_true(result,['representable damping tail ' storage{1} ' ' num2str(count)],rejected,...
                         'a physically representable discarded tail remains an explicit range error');
        if count==256
            result=test_true(result,['collective tail reference ' storage{1}],exp(log(count)+rate)>0,...
                             'individually underflowed exponentials can have a representable collective signal');
        end
    end
end

% Reject subnormal diagonal products while retaining the normal working boundary
for storage={'dense','sparse'}
    for amplitude=[0.5 0.75 1]
        for swapped=[false true]
            initial=amplitude*2^1023*[1;1]; detector=[1;-1];
            if swapped, [initial,detector]=deal(detector,initial); end
            for frequency=[-2^-1074 2^-1074 -2^-1073 2^-1073]
                term=diag([frequency 0]);
                if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(2),{term},initial,detector); rejected=false;
                try
                    pod_evolve(prepared,1,1,1);
                catch caught
                    rejected=contains(caught.message,'normal double working range');
                end
                result=test_true(result,['diagonal product range ' storage{1} ' ' num2str(amplitude) ' ' num2str(swapped) ' ' num2str(frequency)],rejected,...
                                 'a zero or quantised subnormal phase product is rejected before amplification');
            end
            for frequency=[-2^-1021 2^-1021]
                term=diag([frequency 0]);
                if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(2),{term},initial,detector);
                observed=pod_evolve(prepared,1,1,1); expected=-amplitude*pow2(frequency,1023);
                result=test_close(result,['diagonal normal boundary ' storage{1} ' ' num2str(amplitude) ' ' num2str(swapped) ' ' num2str(frequency)],...
                                  imag(observed(2))/expected,1,1e-12,0,...
                                  'a normal component product retains the independently amplified phase');
            end
        end
    end
end

% Check exponent alignment through cancellation, including exact subnormal controls
weak_exp=[-25 -24 -24 -23 -23]; amplitudes=[0.75 0.75 0.5 0.5 0.75];
for storage={'dense','sparse'}
    for swapped=[false true]
        for phase=[1 -1 1i -1i]
            for k=1:5
                initial=phase*[pow2(.5625,512);-pow2(.75,511);pow2(amplitudes(k),weak_exp(k))];
                detector=[pow2(.5,512);pow2(.75,512);pow2(amplitudes(k),weak_exp(k))];
                expected=phase*pow2(amplitudes(k)^2,2*weak_exp(k));
                if swapped, [initial,detector]=deal(detector,initial); expected=conj(expected); end
                term=zeros(3); if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(3),{term},initial,detector);
                if k<=2
                    rejected=false;
                    try
                        pod_evolve(prepared,1,0,0);
                    catch caught
                        rejected=contains(caught.message,'observable product dynamic range');
                    end
                    result=test_true(result,['quantised overlap ' storage{1} ' ' num2str(swapped) ' ' num2str(phase) ' ' num2str(k)],rejected,...
                                     'a nonzero but inexact subnormal alignment cannot certify a finite overlap');
                else
                    observed=pod_evolve(prepared,1,0.3,2);
                    result=test_close(result,['exact subnormal overlap ' storage{1} ' ' num2str(swapped) ' ' num2str(phase) ' ' num2str(k)],...
                                      observed/expected,ones(1,3),1e-12,0,...
                                      'exactly aligned subnormal groups retain the independently known cancellation residual');
                end
            end
        end
    end
end

% Retain independent residuals before filtering equal-exponent product groups
orders=[1 2 3 4;1 3 2 4];
for storage={'dense','sparse'}
    for swapped=[false true]
        for phase=[1 -1 1i -1i]
            for k=1:size(orders,1)
                initial=phase*2^27*[.5;-.5;.875;-.875];
                detector=2^27*[.75;.75-eps(.75);.875;.875]; expected=phase;
                initial=initial(orders(k,:)); detector=detector(orders(k,:));
                if swapped, [initial,detector]=deal(detector,initial); expected=conj(expected); end
                term=zeros(4); if strcmp(storage{1},'sparse'), term=sparse(term); end
                prepared=pod_project(eye(4),{term},initial,detector);
                result=test_close(result,['equal-exponent cancellation ' storage{1} ' ' num2str(swapped) ' ' num2str(phase) ' ' num2str(k)],...
                                  pod_evolve(prepared,1,0.3,2),expected*ones(1,3),1e-12,0,...
                                  'computed dyadic products have an independently exact residual of one before phase rotation');
            end
        end
    end
end

% Preserve projection-slice residuals before long-time phase amplification
orders=[1 2 3 4;1 3 2 4];
for storage={'dense','sparse'}
    for phase=[1 -1 1i -1i]
        for k=1:size(orders,1)
            subspace=.5*ones(4,1); values=[.375,-(.375-eps(.375)),.765625,-.765625];
            term=phase*diag(values(orders(k,:)));
            if strcmp(storage{1},'sparse'), subspace=sparse(subspace); term=sparse(term); end
            prepared=pod_project(subspace,{term},subspace,subspace);
            result=test_close(result,['projection slice residual ' storage{1} ' ' num2str(phase) ' ' num2str(k)],...
                              full(prepared.terms{1})/2^-56,phase,1e-12,0,...
                              'independent dyadic slice contributions retain the projected 2^-56 residual');
            result=test_close(result,['projected long phase ' storage{1} ' ' num2str(phase) ' ' num2str(k)],...
                              pod_evolve(prepared,1,2^56,1),[1 exp(-1i*phase)],1e-12,0,...
                              'a finite long acquisition resolves the independently known nonzero generator');
        end
    end
end

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
         @()pod_basis({raw;NaN(3,2)},[1;0],2,exact_opts),...
         @()pod_basis({raw},1,4,exact_opts),...
         @()pod_project(2*subspace,terms,rho,coil),...
         @()pod_project(subspace,{eye(2)},[1;0;0],[1;0;0]),...
         @()pod_evolve(model,[1 1],0.1,2),...
         @()pod_evolve(bad_model,coefficients,0.1,2),...
         @()pod_evolve(model,coefficients,-0.1,2),...
         @()pod_evolve(model,coefficients,0.1,1.5),...
         @()pod_evolve(model,coefficients,0.1,uint8(255)),...
         @()pod_evolve(model,coefficients,single(0.1),2)};
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

