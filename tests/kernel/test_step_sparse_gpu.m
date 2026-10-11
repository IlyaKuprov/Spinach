% Tests small-Hilbert sparse GPU propagation against spin rotation identities.
% Syntax:
%
%                     result=test_step_sparse_gpu()
%
% Outputs:
%
%     result - CPU controls and GPU sparse/full propagation comparisons
%
% Actual E8/E8 operators reproduce the ideal hard-pulse fixture. Real and
% complex Hermitian generators, promoted and existing GPU input, and matrix
% and cell states are exercised. GPU checks explicitly skip if unavailable.
%
% talos@spindynamics.org

function result=test_step_sparse_gpu()

% Describe the production representation contract
result=new_test_result('kernel/step_sparse_gpu',...
                       'Small Hilbert sparse GPU propagation',...
                       'Sparse GPU storage must not prevent small Hilbert spin rotations.');

% Build the actual two-electron Hilbert basis without interaction assumptions
sys.magnet=0; sys.isotopes={'E8','E8'};
inter.zeeman.scalar={2.002319,2.002319};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
H=operator(spin_system,'Lx','electrons');
K=operator(spin_system,'Ly','electrons');
rho=state(spin_system,'Lz','E8');
assert(issparse(H)&&size(H,1)<spin_system.tols.small_matrix);

% Report unavailable GPU coverage without suppressing production failures
gpu_ok=(exist('canUseGPU','file')~=0)&&canUseGPU();
if ~gpu_ok
    result.messages{end+1}='SKIP: sparse GPU step checks require a usable GPU and Parallel Computing Toolbox.';
end

% Compare both Hermitian generators with independent angular-momentum rotations
for n=1:2
    if n==1
        L=H; angle=pi/2; reference=-full(K);
    else
        L=H+0.7*K; angle=0.37; omega=sqrt(1+0.7^2);
        reference=cos(omega*angle)*full(rho)+...
                  (sin(omega*angle)/omega)*full(0.7*H-K);
    end
    P=expm(-1i*full(L)*angle);
    result=test_close(result,sprintf('CPU analytic fixture=%d',n),...
                      P*full(rho)*P',reference,1e-11,1e-11,...
                      'independent spin rotations fix the sign and magnitude');
    for k=1:(2+4*double(gpu_ok))
        spin_system.sys.enable={};
        if k==3||k==4, spin_system.sys.enable={'gpu'}; end
        if mod(k,2)==1, generator=L; else, generator=full(L); end
        initial=rho;
        if k>=5
            generator=gpuArray(generator); initial=gpuArray(initial);
        end
        for cell_state=0:1
            if cell_state
                state_in={initial,2*initial}; refs={reference,2*reference};
            else
                state_in=initial; refs={reference};
            end
            observed=step(spin_system,generator,state_in,angle);
            label=sprintf('fixture=%d route=%d cell=%d',n,k,cell_state);
            if ~iscell(observed), observed={observed}; end
            for m=1:numel(observed)
                result=test_close(result,[label ' state ' num2str(m)],...
                                  gather(observed{m}),refs{m},1e-11,1e-11,...
                                  'all storage routes retain the same two-sided spin rotation');
            end
        end

        % Preserve input values, storage, and device at zero duration
        zero=step(spin_system,generator,initial,0);
        result=test_true(result,sprintf('fixture=%d route=%d zero',n,k),...
                         isequaln(zero,initial)&&issparse(zero)==issparse(initial)&&...
                         isa(zero,'gpuArray')==isa(initial,'gpuArray'),...
                         'zero-time propagation returns before GPU or storage conversion');
    end
end

end


