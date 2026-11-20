% Checks stimulated-echo ESEEM against dense Hilbert-space propagation
% for an electron and a proton, with and without anisotropic hyperfine
% coupling. Syntax:
%
%                     result=test_eseem_3p()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% Talos <talos@spindynamics.org>

function result=test_eseem_3p()

% Set up the comparison
result=new_test_result('kernel/eseem_3p','Stimulated-echo ESEEM',...
                       'Dense two-spin echo amplitudes and timing boundaries.');
sys.isotopes={'E','1H'};
sys.magnet=2*pi*14.8e6/spin('1H');
bas.approximation={'none'};

% Specify the prepared transverse experiment and a one-orientation grid
parameters.spins={'E'};
parameters.tau=104e-9;
parameters.timestep=4e-9;
parameters.npoints=9;
parameters.grid=struct('alphas',0,'betas',0,'gammas',0,'weights',1);
parameters.serial=true;
receiver=exp(1i*pi/7);

% Express the reference in the ordinary electron-proton product basis
lx=[0 1; 1 0]/2;
lz=diag([1 -1])/2;
Ex=kron(lx,eye(2));
Ez=kron(lz,eye(2));
Nx=kron(eye(2),lx);
Nz=kron(eye(2),lz);
same_e=kron(eye(2),ones(2));

% Compare both types of pulse in both complete Liouville bases
for formalism={'sphten-liouv','zeeman-liouv'}

    % Match the spherical-tensor basis's unit-normalised identity convention
    bas.formalism=formalism{1};
    receiver_scale=1;
    if strcmp(bas.formalism,'sphten-liouv')
        receiver_scale=1/size(Ex,1);
    end

    % Build each spin system from the native tensor input
    for coupling=[0 1]

        % Use the axial tensor tilted by 45 degrees about y
        axis=[1; 0; 1]/sqrt(2);
        inter.coupling.matrix={[],coupling*1e6*(3*(axis*axis')-eye(3)); [],[]};
        spin_system=test_spin_system(sys,inter,bas);
        parameters.rho0=state(spin_system,'Lx','E');
        parameters.coil=receiver*coil_state(spin_system,'Lx','E');
        parameters.pulse_op=operator(spin_system,'Lx','E');
        parameters.offset=coupling*2.75e6;
        H=2*pi*(parameters.offset*Ez-14.8e6*Nz+coupling*1e6*Ez*(1.5*Nx+0.5*Nz));
        U=expm(-1i*H*parameters.tau);

        % Compare raw complex echoes before any processing
        for duration=[0 25e-9]

            % Calculate the native echo
            parameters.pulse_duration=duration;
            fid=powder(spin_system,@eseem_3p,parameters,'esr');

            % Form the pulse and project the reference onto electron-diagonal blocks
            P=expm(-1i*(H*duration+(pi/2)*Ex));
            rho=(P*U*Ex*U'*P').*same_e;
            reference=zeros(parameters.npoints,1);
            for n=1:parameters.npoints

                % Propagate each storage delay independently
                V=U*P*expm(-1i*H*((n-1)*parameters.timestep));
                reference(n)=receiver_scale*conj(receiver)*trace(Ex*V*rho*V');

            end

            % Require agreement in both quadratures and the uncoupled null
            label=sprintf('%s, hyperfine %d, pulse %.0f ns',formalism{1},coupling,duration*1e9);
            result=test_close(result,label,fid,reference,1e-11,0,...
                              'Both quadratures must match the chronological density matrix calculation.');
            if coupling==0

                % On resonance an uncoupled transverse electron is rejected by the filter
                result=test_close(result,[label ' null'],fid,zeros(size(fid)),1e-11,0,...
                                  'An on-resonance transverse electron cannot contribute to electron order zero.');

            end

        end

    end

end

% Check that a single echo samples storage time zero
parameters.npoints=1;
first_echo=powder(spin_system,@eseem_3p,parameters,'esr');
result=test_close(result,'single echo',first_echo,reference(1),1e-11,0,...
                  'Requesting one point must still acquire the zero-storage echo.');

% Check rejection of a negative rectangular pulse duration
parameters.pulse_duration=-25e-9;
try
    powder(spin_system,@eseem_3p,parameters,'esr');
    refused=false;
catch err
    refused=strcmp(err.identifier,'MATLAB:expectedNonnegative');
end
result=test_true(result,'negative pulse duration',refused,...
                 'Pulse duration must be nonnegative.');

end


