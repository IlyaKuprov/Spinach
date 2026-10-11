% Tests GPU norm values and propagation after sparse norm evaluation.
% CPU-only hosts check the CPU norm contract and report GPU coverage absent.
% GPU hosts also exercise the physical state_transfer_wf derivative auxiliaries.
%
% talos@spindynamics.org

function result=test_gpu_cheap_norm()

% Define the norm and propagation contract
result=new_test_result('kernel/gpu_cheap_norm',...
                       'Sparse GPU norm and propagation',...
                       'Norm evaluation must preserve subsequent exponential actions.');
A=sparse([1 -2 3;4 5 -6]);
result=test_close(result,'CPU one-norm',cheap_norm(A),norm(A,1),0,0,...
                  'CPU inputs retain the one-norm contract');
if gpuDeviceCount('available')==0
    result.messages{end+1}='GPU coverage not run: no available device.';
    return;
end

% Check rectangular matrices, vectors, empty shapes, and complex storage
inputs={A,single(A),A+1i*A,sparse([1 -2 3]),sparse([1;-2;3]),...
        sparse(0,3),sparse(3,0),sparse(3,4)};
for n=1:numel(inputs)
    expected=norm(inputs{n},inf);
    actual=gather(cheap_norm(gpuArray(inputs{n})));
    result=test_close(result,['GPU norm ' num2str(n)],actual,expected,1e-13,1e-13,...
                      'GPU inputs retain their matrix or vector infinity-norm');
end

% Preserve native norm handling of non-finite entries
inputs={sparse([NaN 0;1 2]),sparse([Inf 0;1 2])};
for n=1:numel(inputs)
    actual=gather(cheap_norm(gpuArray(inputs{n})));
    result=test_true(result,['GPU non-finite norm ' num2str(n)],...
                     isequaln(actual,norm(inputs{n},inf)),...
                     'NaN and infinity retain the native norm convention');
end

% Build the supplied four-spin wavefunction example's physical generators
sys.magnet=14.1;
sys.isotopes={'1H','1H','13C','13C'};
inter.zeeman.scalar={1.5,2.0,30.0,40.0};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=7.0;
inter.coupling.scalar{1,3}=150;
inter.coupling.scalar{2,4}=150;
inter.coupling.scalar{3,4}=50;
bas.formalism='zeeman-wavef';
bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
H=hamiltonian(assume(spin_system,'nmr'));
parameters.spins={'1H','13C'};
parameters.offset=[1050 5285];
H=frqoffset(spin_system,H,parameters);
controls={operator(spin_system,'Lx','1H'),operator(spin_system,'Ly','1H'),...
          operator(spin_system,'Lx','13C'),operator(spin_system,'Ly','13C')};
rng(101,'twister'); guess=randn(4,100)/10;
for n=1:4, H=H+2*pi*460*guess(n,1)*controls{n}; end
rho=zeros(16,1); rho(end)=1;
zero_mat=sparse(16,16); zero_vec=zeros(16,1); D=controls{1};
generators={H,[H D;zero_mat H],...
            [H D zero_mat;zero_mat H D;zero_mat zero_mat H]};
states={rho,[zero_vec;rho],[zero_vec;zero_vec;rho]};
spin_system.sys.enable={'gpu'};

% Compare the production action with independent dense exponentials
for n=1:numel(generators)
    L=generators{n}; initial=states{n}; dt=1.5e-4;
    expected=expm((-1i*dt)*full(L))*initial;
    actual=gather(step(spin_system,L,initial,dt));
    result=test_close(result,['GPU auxiliary ' num2str(n)],actual,expected,1e-13,1e-13,...
                      'Hermitian and non-Hermitian derivative generators retain the same action');
end

% Retain the already-resident GPU input route for the derivative auxiliary
L=gpuArray(generators{2}); initial=gpuArray(states{2});
expected=expm((-1i*dt)*full(generators{2}))*states{2};
actual=step(spin_system,L,initial,dt);
result=test_true(result,'GPU result residency',isa(actual,'gpuArray'),...
                 'Already-resident operands remain on the GPU');
result=test_close(result,'resident GPU auxiliary',gather(actual),expected,1e-13,1e-13,...
                  'Already-resident non-Hermitian inputs retain the same action');

% Exercise multiple substeps, dissipation, stacked states, and zero duration
rng(1701,'twister');
H=sprandn(64,64,0.2)+1i*sprandn(64,64,0.2); H=(H+H')/2;
L=H-1i*spdiags(linspace(0,1,64)',0,64,64)+spdiags(ones(64,1),1,64,64);
rho=randn(64,3)+1i*randn(64,3);
generators={H,L,L}; states={rho(:,1),rho(:,1),rho};
for n=1:numel(generators)
    expected=expm((-1i*0.3)*full(generators{n}))*states{n};
    actual=gather(step(spin_system,generators{n},states{n},0.3));
    result=test_close(result,['GPU multi-step ' num2str(n)],actual,expected,1e-13,1e-13,...
                      'General sparse generators and stacked states retain their exponential action');
end
result=test_close(result,'zero-time action',step(spin_system,L,rho,0),rho,0,0,...
                  'Zero duration leaves the input unchanged');

end


