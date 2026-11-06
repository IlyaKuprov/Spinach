% Tests Floquet-Magnus averaging against native Fokker-Planck generators.
% Syntax:
%
%                        result=test_fme_context()
%
% Outputs:
%
%     result - native-context integration regression checks
%
% Reuses rotor_stack() for Fourier input and drifts() with singlerot()
% for an independent two-spin MAS reference. The enlarged generator's
% spatial derivative is never included in the Fourier input. The native
% positive-rate transport has physical angular frequency -2*pi*rate.
%
% talos@spindynamics.org

function result=test_fme_context()

% Describe the native-context integration target
result=new_test_result('kernel/fme_context','Floquet-Magnus native FP reference',...
                       'averaging must match the existing drifts and singlerot conventions');

% Build a noncommuting two-spin dipolar system in full Liouville space
sys.magnet=9.4; sys.isotopes={'1H','1H'};
inter.zeeman.scalar={5,-2};
inter.coordinates={[0 0 0],[0.6 0.8 1.4]};
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
rho=state(spin_system,'Lx','all'); rho=rho/norm(rho);
coil=state(spin_system,'L+','all'); coil=coil/norm(coil);

% Use the native context to generate both the rotor transport and interactions
parameters.spins={'1H'}; parameters.offset=0;
parameters.rate=5e5; parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='single_crystal'; parameters.serial=true;
parameters.orientation=[0 0 0]; parameters.masframe='rotor';
parameters.rframes={};
parameters.verbose=0; omega=-2*pi*parameters.rate;
duration=0.73/parameters.rate; dim=numel(rho); signals=zeros(1,2);
for k=1:2

    % Compare two independently refined native spatial Fourier cutoffs
    parameters.max_rank=6+2*k;
    [ensemble,spc_dim]=drifts(spin_system,@singlerot,parameters,'nmr');
    L=ensemble{1}{1};
    [stack,phases]=rotor_stack(spin_system,parameters,'nmr');
    initial=kron(sparse(1,1,1,spc_dim,1),rho);
    detector=kron(ones(spc_dim,1),coil);
    signals(k)=detector'*expm(-1i*full(L)*duration)*initial;

    % Obtain the five harmonics from native spin-space rotor Hamiltonians
    H=repmat({sparse(dim,dim)},1,5);
    for n=-2:2
        for q=1:spc_dim
            H{n+3}=H{n+3}+stack{q}*exp(-1i*n*phases(q))/spc_dim;
        end
    end

    % Check the block reconstruction before interpreting a numerical comparison
    for q=1:spc_dim
        reconstructed=sparse(dim,dim);
        for n=-2:2
            reconstructed=reconstructed+H{n+3}*exp(1i*n*phases(q));
        end
        idx=(q-1)*dim+(1:dim);
        result=test_close(result,'native harmonic reconstruction',reconstructed,...
                          L(idx,idx),1e-12,1e-9,...
                          'a secular dipolar MAS Hamiltonian contains only ranks zero to two');
    end

    % Propagate the initial zero-phase slice with both endpoint kicks
    [F,K]=floquet_magnus(H,omega,3,0,'van_vleck');
    rank=(numel(K)-1)/2; start=sparse(dim,dim); finish=sparse(dim,dim);
    for n=-rank:rank
        start=start+K{n+rank+1};
        finish=finish+K{n+rank+1}*exp(1i*n*omega*duration);
    end
    evolved=expm(-1i*full(finish))*expm(-1i*full(F)*duration)*...
            expm(1i*full(start))*rho;
    signal=coil'*evolved;
    result=test_close(result,'native FP complex signal',signal,signals(k),...
                      1e-9,1e-9,'the complete kick approximation must follow native signed rotor transport');

end

% The spatial discretisation must be converged independently of averaging order
result=test_close(result,'native FP rank convergence',signals(1),signals(2),...
                  1e-9,1e-9,'increasing the FP rank must not change the reference signal');

end


