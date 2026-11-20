% Checks repeated DNP transfer against chronological propagation of a
% detuned electron-proton density matrix. Syntax:
%
%                      result=test_dnp_train()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% Talos <talos@spindynamics.org>

function result=test_dnp_train()

% Set up a small complete Liouville basis
sys.magnet=2*pi*14.8e6/spin('1H');
sys.isotopes={'E','1H'};
bas.formalism='zeeman-liouv';
bas.approximation={'none'};
spin_system=test_spin_system(sys,struct(),bas);

% Define a detuned drift with a pseudosecular hyperfine term
sx=[0 1; 1 0]/2;
sy=[0 -1i; 1i 0]/2;
sz=diag([1 -1])/2;
Ex=kron(sx,eye(2));
Ey=kron(sy,eye(2));
Iz=kron(eye(2),sz);
H0=2*pi*(8e6*kron(sz,eye(2))-14.8e6*Iz+...
         1.5e6*kron(sz,sx)+0.5e6*kron(sz,sz));
H=kron(eye(4),H0)-kron(transpose(H0),eye(4));
R=-2e4*speye(16);
K=-3e4*speye(16);

% Specify unequal slices with noncommuting controls and a complex receiver
receiver=-Iz+1i*Ey;
parameters.rho0=state(spin_system,'Lx','E');
parameters.coil=receiver(:);
parameters.pulse_ops={operator(spin_system,'Lx','E'),...
                      operator(spin_system,'Ly','E')};
parameters.slice_durs=[5 9 7]*1e-9;
parameters.amplitudes={2*pi*[12 0 -16]*1e6,2*pi*[0 20 8]*1e6};

% Build independent reference signals by reapplying each physical slice
reference=zeros(8,1);
rho=Ex;
reference(1)=trace(receiver'*rho);
for n=1:7

    % Evolve in Hilbert space, including the scalar Liouville loss
    for k=1:numel(parameters.slice_durs)

        % Apply the current segment on both sides of the density matrix
        dt=parameters.slice_durs(k);
        Hslice=H0+parameters.amplitudes{1}(k)*Ex+parameters.amplitudes{2}(k)*Ey;
        U=expm(-1i*dt*Hslice);
        rho=exp(-5e4*dt)*(U*rho*U');

    end

    % Read the signed complex transfer without normalising
    reference(n+1)=trace(receiver'*rho);

end

% Compare one period before checking two and seven periods
result=new_test_result('kernel/dnp_train','Repeated DNP train',...
                       'Checks chronological slice order and complex detection across repeats.');
for n=[1 2 7]

    % Compare all period boundaries, including the unirradiated point
    parameters.nloops=n;
    signal=dnp_train(spin_system,parameters,H,R,K);
    result=test_close(result,sprintf('%d period(s)',n),signal,reference(1:n+1),1e-10,0,...
                      'The train must match direct complex detection of the propagated density matrix.');

end

% Reject a waveform with a missing Cartesian amplitude
parameters.amplitudes{2}(end)=[];
try

    % Call the production entry point with inconsistent input
    dnp_train(spin_system,parameters,H,R,K);
    refused=false;

catch err

    % Check the input-specific rejection
    refused=strcmp(err.identifier,'Spinach:dnp_train:badAmplitudes');

end
result=test_true(result,'inconsistent slice counts',refused,...
                 'Each Cartesian channel needs exactly one amplitude per slice.');

end


