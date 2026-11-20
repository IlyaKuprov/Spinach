% Tests shifted-echo STMAS time ordering, complex detection, and
% frequency calibration against explicit spin-3/2 matrices.
% Syntax:
%
%                       result=test_stmas()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_stmas()

% Define the test and construct the native basis
result=new_test_result('kernel/stmas','Shifted-echo STMAS',...
                       'the complex echo agrees with independently ordered finite pulses.');
sys.magnet=0;
sys.isotopes={'87Rb'};
sys.parallel={'processes',1};
sys.output='hush';
inter.zeeman.scalar={0};
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);

% Specify a small complex echo with two nonzero indirect increments
parameters.rate=50e3;
parameters.spc_dim=1;
parameters.npoints=[2 8];
parameters.pulse_dur=[1.6 1.6 25]*1e-6;
parameters.pulse_amp=2*pi*[100e3 100e3 10e3];
parameters.tau=80e-6;
parameters.rho0=state(spin_system,'Lz',1);
parameters.coil=coil_state(spin_system,'CT-',1);
period=20e-6;
sat_times=[18.4 38.4]*1e-6;
ct_times=[96.3555555555556 114.1333333333333]*1e-6;
direct_times=(0:7)*period;

% Set independent Hilbert operators, masks, and receiver
lz=diag([1.5 0.5 -0.5 -1.5]);
lp=diag([sqrt(3) 2 sqrt(3)],1);
lx=(lp+lp')/2;
receiver=zeros(4);
receiver(3,2)=1;
sat_support=logical(diag([1 0 1],1));
ct_support=logical(receiver');

% Compare a resolved quadrupole and its zero-interaction limit
for quad_hz=[2100 0]

    % Build the generator and run the native pulse callback
    drift=2*pi*(4700*lz+quad_hz*(lz^2-(5/4)*eye(4)));
    H=hilb2liouv(drift,'comm');
    fid=stmas(spin_system,parameters,H,sparse(16,16),sparse(16,16));

    % Build finite pulse propagators in Hilbert space
    pulse=cell(1,3);
    for n=1:3
        pulse{n}=expm(-1i*(drift+parameters.pulse_amp(n)*lx)*parameters.pulse_dur(n));
    end

    % Calculate the two complex traces using entry selection in Hilbert space
    expected=complex(zeros(8,2));
    for n=1:2

        % Propagate the selected satellites and central coherence
        rho=(pulse{1}*lz*pulse{1}').*sat_support;
        free=expm(-1i*drift*sat_times(n));
        rho=(pulse{2}*free*rho*free'*pulse{2}').*ct_support;
        free=expm(-1i*drift*ct_times(n));
        rho=pulse{3}*free*rho*free'*pulse{3}';

        % Reference detection, including the sample immediately after inversion
        for k=1:8
            free=expm(-1i*drift*direct_times(k));
            expected(k,n)=trace(receiver'*free*rho*free');
        end

    end

    % Check the time-domain echo and both transformed axes
    label=[num2str(quad_hz) ' Hz quadrupole'];
    result=test_true(result,[label ' finite signal'],norm(expected,'fro')>0.01,...
                     'the finite pulse reference must distinguish phase and amplitude');
    result=test_close(result,[label ' complex echo'],fid,expected,1e-7,1e-7,...
                      'the two increments follow signed CT detection and the shifted echo clock');
    result=test_close(result,[label ' complex spectrum'],fftshift(fft2(fid)),...
                      fftshift(fft2(expected)),4e-7,1e-7,...
                      'Fourier processing preserves direct and indirect dimensions and complex phase');

end

% Refuse a negative first free-delay interval
parameters.pulse_dur=[30 20 25]*1e-6;
try
    stmas(spin_system,parameters,H,sparse(16,16),sparse(16,16));
    refused=false;
catch err
    refused=contains(err.message,'negative first satellite free delay');
end
result=test_true(result,'impossible pulse timing',refused,...
                 'the first indirect increment cannot overlap excitation and conversion');

end


