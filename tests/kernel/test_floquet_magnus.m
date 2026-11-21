% Tests finite-harmonic Floquet-Magnus averaging and periodic kicks. Syntax:
%
%                       result=test_floquet_magnus()
%
% Outputs:
%
%     result - regression checks, failures, and explanatory messages
%
% Covers exact commuting limits, time-origin covariance, degenerate
% eigenspaces, input boundaries, and direct time-domain convergence of
% both gauges through inverse-frequency power three. Numerical reference
% propagation uses independently converged ode113 integration.
%
% talos@spindynamics.org

function result=test_floquet_magnus()

% Create a result using the standard regression harness
result=new_test_result('kernel/floquet_magnus','Harmonic Floquet-Magnus expansion',...
                       'effective generators and kicks must reproduce periodic dynamics');

% Check an exactly degenerate constant generator in both gauges
H={sparse(diag([1 1 3]))};
for theory={'van_vleck','aht'}
    for order=0:3
        [F,K]=floquet_magnus(H,-7,order,0.13,theory{1});
        result=test_close(result,'constant generator',F,H{1},0,0,...
                          'a constant generator has no averaging corrections');
        result=test_close(result,'constant kick',K{1},zeros(3),0,0,...
                          'a constant generator has no micromotion');
    end
end

% Check a commuting waveform and its degenerate eigenspace projector
rotation=expm(1i*[0 1 2;1 0 3;2 3 0]);
A=rotation*diag([1 1 3])*rotation';
H={0.2i*A,0.3*A,A,0.3*A,-0.2i*A};
omega=17; t0=0.021; tend=t0+0.63*2*pi/omega;
[F,K]=floquet_magnus(H,omega,3,t0,'van_vleck');
phase=(tend-t0)+0.6*(sin(omega*tend)-sin(omega*t0))/omega-...
      0.2*(cos(2*omega*tend)-cos(2*omega*t0))/omega;
P=kick_propagator(F,K,omega,t0,tend);
result=test_close(result,'commuting exact propagation',P,expm(-1i*A*phase),...
                  2e-13,2e-13,'commuting Fourier components integrate exactly');
[vectors,values]=eig(full(F));
[projector,~]=qr(vectors(:,abs(diag(values)-1)<1e-10),0);
projector=projector*projector';
result=test_close(result,'degenerate spectral projector',projector,...
                  rotation(:,1:2)*rotation(:,1:2)',2e-12,2e-12,...
                  'a degenerate eigenspace is compared without choosing eigenvectors');

% Preserve the caller random stream while constructing a reproducible family
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng));
rng(41); H=cell(1,5);
H{3}=randn(3)+1i*randn(3); H{3}=(H{3}+H{3}')/2;
for n=4:5
    H{n}=(randn(3)+1i*randn(3))/3; H{6-n}=H{n}';
end

% Check order-dependent errors against an independent non-stroboscopic reference
for theory={'van_vleck','aht'}
    errors=zeros(4,2);
    for k=1:2
        omega=20*2^(k-1); period=2*pi/omega;
        t0=0.17*period; tend=t0+0.73*period;
        reference=direct_propagator(H,omega,t0,tend,1e-12);
        refined=direct_propagator(H,omega,t0,tend,2e-13);
        result=test_close(result,'reference convergence',reference,refined,...
                          2e-11,2e-11,'tightening the ODE tolerance must not change the reference');
        for order=0:3
            [F,K]=floquet_magnus(H,omega,order,t0,theory{1});
            P=kick_propagator(F,K,omega,t0,tend);
            errors(order+1,k)=norm(P-refined,'fro');
            result=test_close(result,'Hermitian effective generator',F,F',...
                              1e-12,1e-12,'Hermitian Fourier input must remain Hermitian');
            kick=kick_value(K,omega,tend);
            result=test_close(result,'Hermitian kick',kick,kick',1e-12,1e-12,...
                              'the physical kick is Hermitian without post-processing');
            result=test_close(result,'unitary propagator',P'*P,eye(3),...
                              1e-12,1e-12,'Hermitian input gives unitary factorised propagation');
            if strcmp(theory{1},'aht')
                kick=kick_value(K,omega,t0);
                result=test_close(result,'AHT time origin',kick,zeros(3),...
                                  1e-12,1e-12,'the AHT kick vanishes at its specified time origin');
            else
                result=test_close(result,'Van Vleck zero mean',K{(numel(K)+1)/2},...
                                  zeros(3),0,0,'the Van Vleck kick has zero Fourier mean');
            end
        end
    end
    for order=0:3
        slope=log2(errors(order+1,1)/errors(order+1,2));
        result=test_true(result,sprintf('%s order %d scaling',theory{1},order),...
                         slope>order+1.7,...
                         sprintf('measured power %.3f must approach %d',slope,order+2));
    end
end

% Check Fourier-phase covariance without confusing it with gauge equivalence
omega=31; shift=0.037; shifted=H;
for n=-2:2
    shifted{n+3}=H{n+3}*exp(1i*n*omega*shift);
end
for theory={'van_vleck','aht'}
    [F,K]=floquet_magnus(H,omega,3,shift,theory{1});
    [shifted_f,shifted_k]=floquet_magnus(shifted,omega,3,0,theory{1});
    result=test_close(result,'time-origin effective generator',shifted_f,F,...
                      1e-12,1e-12,'a shift of Fourier phases and time origin describes the same drive');
    for n=1:numel(K)
        harmonic=n-(numel(K)+1)/2;
        result=test_close(result,'time-origin kick coefficient',shifted_k{n},...
                          K{n}*exp(1i*harmonic*omega*shift),1e-12,1e-12,...
                          'kick coefficients transform with the same Fourier phases');
    end
end

% Test stroboscopic evolution and sensitivity to removing the kick
omega=40; t0=0.17*2*pi/omega; tend=t0+2*pi/omega;
[F,K]=floquet_magnus(H,omega,3,t0,'van_vleck');
reference=direct_propagator(H,omega,t0,tend,2e-13);
complete=kick_propagator(F,K,omega,t0,tend);
bare=expm(-1i*full(F)*(tend-t0));
result=test_true(result,'kick is necessary',norm(complete-reference,'fro')<...
                 norm(bare-reference,'fro')/10,...
                 'omitting micromotion must degrade this noncommuting benchmark');
result=test_true(result,'complex observable detects missing kick',...
                 abs(complete(1,2)-reference(1,2))<abs(bare(1,2)-reference(1,2)),...
                 'eigenvalues alone cannot detect an incorrect eigenframe');

% Check negative angular frequency against direct propagation
omega=-40; t0=0.013; tend=t0+0.7*2*pi/abs(omega);
[F,K]=floquet_magnus(H,omega,3,t0,'van_vleck');
reference=direct_propagator(H,omega,t0,tend,2e-13);
P=kick_propagator(F,K,omega,t0,tend);
result=test_true(result,'negative frequency',norm(P-reference,'fro')<2e-4,...
                 'the signed frequency controls both denominators and Fourier phases');

% Reject malformed caller inputs without silently adapting their shapes
bad_calls={@()floquet_magnus(H.',40,1,0,'aht'),...
           @()floquet_magnus(H,0,1,0,'aht'),...
           @()floquet_magnus(H,40,4,0,'aht'),...
           @()floquet_magnus(H,40,1,Inf,'aht'),...
           @()floquet_magnus(H,40,1,0,'unknown'),...
           @()floquet_magnus(H,int32(40),1,0,'aht'),...
           @()floquet_magnus(H,40,int32(1),0,'aht'),...
           @()floquet_magnus(H,40,1,int32(0),'aht'),...
           @()floquet_magnus({eye(2),eye(3),eye(2)},40,1,0,'aht')};
for n=1:numel(bad_calls)
    caught=false;
    try
        bad_calls{n}();
    catch
        caught=true;
    end
    result=test_true(result,'input rejection',caught,...
                     'invalid expansion controls or coefficient shapes must be rejected');
end

% Restore the caller stream before returning the test result
clear restore_rng

end

% Evaluate a finite Fourier series for the periodic kick
function value=kick_value(K,omega,time)

% Sum every supplied harmonic without FFT aliasing
rank=(numel(K)-1)/2; value=zeros(size(K{1}));
for n=-rank:rank
    value=value+K{n+rank+1}*exp(1i*n*omega*time);
end

end

% Form the complete effective-Hamiltonian and kick approximation
function P=kick_propagator(F,K,omega,t0,tend)

% Keep both endpoint kicks in their defining order
start=kick_value(K,omega,t0); finish=kick_value(K,omega,tend);
P=expm(-1i*full(finish))*expm(-1i*full(F)*(tend-t0))*expm(1i*full(start));

end

% Integrate the original Schrodinger equation independently of the recurrence
function P=direct_propagator(H,omega,t0,tend,tolerance)

% Use controlled ODE tolerances rather than an effective generator
opts=odeset('RelTol',tolerance,'AbsTol',tolerance/10);
dim=size(H{1},1);
[~,trajectory]=ode113(@(time,state)ode_rhs(time,state,H,omega),...
                     [t0 tend],reshape(eye(dim),dim^2,1),opts);
P=reshape(trajectory(end,:),dim,dim);

end

% Compute the original time-dependent matrix action
function derivative=ode_rhs(time,state,H,omega)

% Construct the Hamiltonian from the unmodified input harmonics
rank=(numel(H)-1)/2; dim=size(H{1},1); ham=zeros(dim);
for n=-rank:rank
    ham=ham+H{n+rank+1}*exp(1i*n*omega*time);
end
derivative=reshape(-1i*ham*reshape(state,dim,dim),dim^2,1);

end


