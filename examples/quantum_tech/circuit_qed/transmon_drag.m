% DRAG correction of a resonant Gaussian pulse on a three-level
% Duffing transmon in the laboratory frame. The derivative of the
% envelope, scaled by the DRAG detuning parameter, is fed into the
% quadrature channel of the IQ mixer; this suppresses the leakage
% into the second excited state and improves the fidelity of the
% 90 degree rotation on the qubit subspace. A further improvement
% comes from numerical optimisation of the mixer parameters. The
% rows of pulse_params are the plain Gaussian, the analytically
% corrected DRAG, and the numerically optimised DRAG pulse; the
% columns are the envelope amplitude, the DRAG detuning, the lo-
% cal oscillator frequency, the mixer phase, and the DRAG quad-
% rature switch. The propagator is taken into the frame of the
% drift Hamiltonian before it is compared with the target gate,
% which is specified in that frame; without this the comparison
% would only be meaningful when every transmon level phase hap-
% pens to close on a multiple of two pi over the pulse duration.
% Model and parameters from the DRAG example of the paraqeet
% package.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_drag()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Transmon frequency and anharmonicity
frq_01=4.8e9; anharm=-200e6;
inter.modes.frqs={frq_01};
inter.modes.anharms={anharm};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Laboratory frame drift Hamiltonian
H0=hamiltonian(assume(spin_system,'labframe'));

% Transmon quadrature drive operator
Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
Dx=full(Cr+An); H0=full(H0);

% Pulse duration and Gaussian envelope width
t_dur=20e-9; sigma=t_dur/8;

% Analytical DRAG detuning parameter for a Duffing transmon
drag_det=2*(2*pi)*anharm;

% Rows: plain Gaussian, analytical DRAG, optimised DRAG
pulse_params=[3.000000000000000e8 drag_det 2*pi*frq_01 0 0;
              3.000000000000000e8 drag_det 2*pi*frq_01 0 1;
              2.455066180403117e8 -2.492899680375753e9 3.015929426382811e10 -1.886372179487061e-4 1];

% Target gate is a 90 degree rotation on the qubit subspace
target=[cos(pi/4) -1i*sin(pi/4) 0; -1i*sin(pi/4) cos(pi/4) 0; 0 0 1];

% Internal time stepping of the source calculation
nsub=2000; dts=1e-11;

% Preallocate fidelities and population trajectories
fids=zeros(1,3); pops=zeros(3,nsub+1,3);

% Loop over the three pulses
for m=1:3

    % Extract mixer parameters
    amp=pulse_params(m,1); delta=pulse_params(m,2);
    w_lo=pulse_params(m,3); phi=pulse_params(m,4);
    drag=pulse_params(m,5);

    % Propagate on the internal midpoint grid
    U=eye(3); psi=[1; 0; 0]; pops(:,1,m)=abs(psi).^2;
    for k=1:nsub

        % Gaussian envelope and its derivative at the midpoint
        tmid=(k-0.5)*dts;
        env=amp*exp(-(tmid-t_dur/2)^2/(2*sigma^2));
        denv=-env*(tmid-t_dur/2)/sigma^2;

        % IQ mixer signal with the DRAG quadrature
        sig=env*cos(w_lo*tmid-phi)-drag*(denv/delta)*sin(w_lo*tmid-phi);

        % Propagate through the slice
        U=expm(-1i*(H0+sig*Dx)*dts)*U;
        pops(:,k+1,m)=abs(U(:,1)).^2;

    end

    % Move the propagator into the frame of the drift Hamiltonian
    U=expm(1i*H0*t_dur)*U;

    % Process fidelity of the gate over all three levels
    fids(m)=abs(mean(diag(target'*U)))^2;

end

% Validate the fidelities against the source calculation
if max(abs(fids-[0.969588 0.973797 0.994239]))>1e-3
    error('pulse fidelities deviate from the reference values.');
end

% Validate the leakage suppression by the DRAG quadrature
if ~(pops(3,end,2)<pops(3,end,1)/2)
    error('DRAG quadrature did not suppress the leakage.');
end

% Plot ground state pulse response with and without DRAG
time_axis=1e9*dts*(0:nsub);
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,pops(:,:,1)','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns'); kylabel('level populations');
ktitle('plain Gaussian pulse');
klegend({'$|0\rangle$','$|1\rangle$','$|2\rangle$'},'Location','Best');
subplot(1,2,2); plot(time_axis,pops(:,:,2)','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns'); kylabel('level populations');
ktitle('analytical DRAG pulse');
klegend({'$|0\rangle$','$|1\rangle$','$|2\rangle$'},'Location','Best');

end

