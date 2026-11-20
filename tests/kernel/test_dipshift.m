% Tests DIPSHIFT echo timing, simultaneous RF, and the matched
% zero-width limit with a separate dense two-spin calculation.
%
% Syntax: result=test_dipshift()
%
% Outputs:
%
%    result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_dipshift()

% Set up a coupled pair with offsets to make finite RF phase-sensitive
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'1H','13C'};
inter.coupling.scalar={[] 7800; [] []};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);

% Specify the rotor, shifts, and complex detection
parameters.spins=sys.isotopes;
parameters.offset=[1100 -2300];
parameters.rate=6000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=2;
parameters.grid='single_crystal';
parameters.serial=true;
parameters.time_shift=[0 0.13 0.37 0.5]/parameters.rate;
parameters.rho0=state(spin_system,'Lx','13C');
parameters.coil=coil_state(spin_system,'L+','13C');

% Compare all three schedules to dense Hilbert-space propagation
result=new_test_result('kernel/dipshift','DIPSHIFT echoes',...
                       'Pulse timing and complex detection under singlerot.');
curves=complex(zeros(3,numel(parameters.time_shift)));
widths=[0 5e-6 0];
labels={'X-only','finite X,Y','zero-width X,Y'};
for n=1:numel(widths)

    % Compute the echo through the production context
    parameters.pulse_width=widths(n);
    parameters.xy_block=(n~=1);
    curves(n,:)=singlerot(spin_system,@dipshift,parameters,'nmr');

    % Use a direct four-level calculation with the literal delay arms
    overlap=parameters.coil'*parameters.rho0;
    ref=dense_echo(parameters,7800);
    result=test_close(result,labels{n},curves(n,:),overlap*ref,1e-9,1e-9,...
                      'DIPSHIFT must preserve the signed complex echo.');

end

% Shrink the X,Y pulse width at fixed flip angle
parameters.pulse_width=5e-10;
short_curve=singlerot(spin_system,@dipshift,parameters,'nmr');
result=test_close(result,'matched zero-width limit',short_curve,curves(3,:),...
                  1e-9,1e-4,'Short X,Y pulses must approach the X,Y ideal schedule.');
result=test_true(result,'distinct hard schedules',...
                 norm(curves(1,:)-curves(3,:))>0.01,...
                 'X-only is not the zero-width reference for the X,Y train.');

% Check the dipole-free, on-resonance negative echo
inter.coupling.scalar={[] []; [] []};
spin_system=test_spin_system(sys,inter,bas);
parameters.offset=[0 0];
parameters.pulse_width=5e-6;
null_curve=singlerot(spin_system,@dipshift,parameters,'nmr');
result=test_close(result,'uncoupled echo',null_curve,...
                  -overlap*ones(size(null_curve)),1e-9,1e-9,...
                  'An observed-channel Y pi pulse sends Lx to minus Lx.');

end

% Evaluate separate arm delays using spin-half matrices and dense exponentials
function signal=dense_echo(parameters,coupling)

% Construct the carrier-offset and weak-coupling Hamiltonian
X=[0 1; 1 0]/2;
Y=[0 -1i; 1i 0]/2;
Z=diag([1 -1])/2;
Iz=kron(Z,eye(2));
Sz=kron(eye(2),Z);
H=2*pi*(coupling*Iz*Sz+parameters.offset(1)*Iz+parameters.offset(2)*Sz);
rho0=kron(eye(2),X);
coil=kron(eye(2),X+1i*Y);
Ix=kron(X,eye(2));
Iy=kron(Y,eye(2));
Sy=kron(eye(2),Y);

% Use the pulse width to form the integrated pulse Hamiltonians
period=1/parameters.rate;
width=parameters.pulse_width;
quarter=period/4;
half=period/2;
exp_drift=@(delay)expm(-1i*H*delay);
pulse_x=expm(-1i*(pi*Ix+H*width));
pulse_y=expm(-1i*(pi*Iy+H*width));
pulse_echo=expm(-1i*(pi*(Ix+Sy)+H*width));
echo_delay=exp_drift(period-width/2);
echo=echo_delay*pulse_echo*echo_delay;

% Apply the two different arm schedules chronologically from the right
signal=complex(zeros(size(parameters.time_shift)));
for n=1:numel(signal)

    % Construct the arms from their separate delays and pulses
    shift=parameters.time_shift(n);
    if parameters.xy_block
        quarter_delay=exp_drift(quarter-width/2);
        block=quarter_delay*pulse_y*quarter_delay^2*pulse_x*quarter_delay;
        out_delay=exp_drift(shift+quarter);
        in_delay=exp_drift(half-shift+quarter);
        first_arm=in_delay*block*out_delay;
        second_arm=out_delay*block*in_delay;
    else
        out_delay=exp_drift(shift);
        in_delay=exp_drift(half-shift);
        half_delay=exp_drift(half);
        first_arm=in_delay*pulse_x*half_delay*out_delay;
        second_arm=out_delay*half_delay*pulse_x*in_delay;
    end

    % Detect after the complete density-matrix evolution
    U=second_arm*echo*first_arm;
    signal(n)=trace(coil'*U*rho0*U')/trace(coil'*rho0);

end

end


