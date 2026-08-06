% Ramsey chevron of a three-level transmon in the Duffing ap-
% proximation. A nominal pi/2 pulse prepares a coherence, and
% detuning during free evolution produces Ramsey fringes.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_ramsey_chevron()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Transmon in the rotating frame
inter.modes.frqs={0};
inter.modes.anharms={-260e6};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Anharmonicity part from the declared interactions
spin_system=assume(spin_system,'cavity');
H0=hamiltonian(spin_system);

% Transmon operators
Cr=operator(spin_system,'C',1);
An=operator(spin_system,'A',1);
N=operator(spin_system,'N',1);

% Free-evolution parameters
detunings=2*pi*linspace(-20e6,20e6,256);
time_axis=linspace(0,1.0e-6,256);

% Nominal pi/2 pulse propagator
X=(Cr+An)/2; U90=expm(-1i*(pi/2)*full(X));

% Initial density matrix
rho=state(spin_system,'BL1',1);
rho=U90*rho*U90';

% Detection state
L2=state(spin_system,'BL2',1);

% Loop over offsets
answer=zeros(numel(detunings),...
             numel(time_axis));
for n=1:numel(detunings)

    % Build the free-evolution Hamiltonian
    H=detunings(n)*N+H0; H=(H+H')/2;

    % Propagate the Ramsey interferometer
    for k=1:numel(time_axis)
        U=expm(-1i*full(H)*time_axis(k));
        rho_t=U90*U*rho*U'*U90';
        answer(n,k)=real(trace(L2*rho_t));
    end

end

% Plot the Ramsey chevron
kfigure(); imagesc(1e6*time_axis,detunings/(2*pi*1e6),answer);
axis xy tight; kxlabel('time, $\mu$s');
kylabel('detuning, MHz'); kcolourbar;
ktitle('transmon Ramsey chevron');

end

