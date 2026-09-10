% Adiabatic preparation of the Stark-dressed states of a transmon
% coupled dispersively to a cavity, following Appendix F of Yunwei
% Lu's PhD thesis (Northwestern University, 2026). The transmon is
% a Duffing oscillator in the frame rotating with the drive; the
% cavity is in its own rotating frame and is coupled to the trans-
% mon through the dispersive shift. The in-phase drive is ramped
% with a Gaussian leading edge from zero to the dynamical sweet-
% spot amplitude of Eq. (4.22), and the DRAG-type quadrature of
% Eqs. (F.14)-(F.16) cancels the leading nonadiabatic coupling.
% The state preparation infidelity of Eq. (F.18) is computed for
% the bare |g,0> and |g,1> states as a function of the ramp time
% for several drive detunings, reproducing the closed-system part
% of Fig. F.1; the thesis figure also includes decoherence, which
% makes its curves turn up at long ramp times.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_stark_ramp()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T4','C3'};

% Transmon anharmonicity, coupling, and transmon-cavity detuning
anharm=-200e6; coupling=100e6; delta_bc=1e9;

% Dispersive shift in the transmon-cavity system
chi=2*(coupling/delta_bc)^2*anharm;

% Transmon-drive detunings and DRAG quadrature switches
delta_bd=[-6e6 -12e6 -19e6 -19e6]; drag=[1 1 1 0];

% Ramp times and the propagation time step
ramp_times=1e-9*(50:10:200); dt=0.5e-9;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Preallocate the infidelities
infid=zeros(numel(delta_bd),numel(ramp_times));

% Loop over the detunings
for n=1:numel(delta_bd)

    % Transmon in the drive frame, cavity in its own frame
    inter.modes.frqs={delta_bd(n) 0};
    inter.modes.anharms={anharm []};

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Rotating frame Hamiltonian with the dispersive coupling
    H0=hamiltonian(assume(spin_system,'labframe'));
    H0=H0+2*pi*chi*operator(spin_system,{'N','N'},{1,2});

    % Transmon drive operators for the two quadratures
    Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
    Dx=Cr+An; Dy=-1i*(An-Cr);

    % Bare |g,0> and |g,1> states from the projector diagonals
    [~,idx_g0]=max(diag(state(spin_system,{'BL1','BL1'},{1,2})));
    [~,idx_g1]=max(diag(state(spin_system,{'BL1','BL2'},{1,2})));
    psi0=full(eye(size(H0,1))); psi0=psi0(:,[idx_g0 idx_g1]);

    % Sweet-spot drive amplitude for the current detuning
    omega_0=sqrt(delta_bd(n)^3/(4*anharm));

    % Dressed eigenstates connected to the bare states at full drive
    [V,~]=eig(full(H0+2*pi*omega_0*Dx)); [~,idx]=max(abs(V'*psi0));
    targets=V(:,idx);

    % Loop over the ramp times
    for k=1:numel(ramp_times)

        % Gaussian edge parameters and the time grid midpoints
        tau=ramp_times(k); sigma=tau/4; nsteps=round(tau/dt);
        tmid=(0.5:1:nsteps)*(tau/nsteps); psi=psi0;

        % In-phase envelope, its derivative, and the DRAG quadrature
        env_shape=exp(-(tmid-tau).^2/(2*sigma^2)); baseline=exp(-tau^2/(2*sigma^2));
        omega_i=omega_0*(env_shape-baseline)/(1-baseline);
        omega_dot=-omega_0*(tmid-tau).*env_shape/(sigma^2*(1-baseline));
        omega_q=-drag(n)*omega_dot/(2*pi*delta_bd(n));

        % Propagate both bare states through the ramp
        for m=1:nsteps
            H=H0+2*pi*omega_i(m)*Dx+2*pi*omega_q(m)*Dy;
            psi=propagator(spin_system,H,tau/nsteps)*psi;
        end

        % Infidelity from the overlaps with the dressed target states
        infid(n,k)=1-mean(abs(diag(targets'*psi)).^2);

    end

end

% Report the infidelities at selected ramp times
for k=[1 6 16]
    fprintf('tau=%3.0f ns: I=%.3e (6 MHz, DRAG), %.3e (12 MHz, DRAG), %.3e (19 MHz, DRAG), %.3e (19 MHz, Gaussian)\n',...
            1e9*ramp_times(k),infid(:,k));
end

% Validate the leakage suppression by the DRAG quadrature
if ~(infid(3,1)<infid(4,1)/10)
    error('DRAG quadrature did not suppress the nonadiabatic leakage.');
end

% Validate the improvement of adiabaticity with the ramp time
if any(infid(:,end)>=infid(:,1))
    error('infidelity does not decrease with the ramp time.');
end

% Plot the infidelities against the ramp time
kfigure(); semilogy(1e9*ramp_times,infid','LineWidth',1.5);
axis tight; kgrid; kxlabel('ramp time, ns'); kylabel('infidelity');
ktitle('Stark-dressed state preparation');
klegend({'$-\Delta_{bd}/2\pi=6$ MHz (DRAG)','$-\Delta_{bd}/2\pi=12$ MHz (DRAG)',...
         '$-\Delta_{bd}/2\pi=19$ MHz (DRAG)','$-\Delta_{bd}/2\pi=19$ MHz (Gaussian)'},...
         'Location','Best');

end

