% Pulsed-field magnetisation of a dimer of two S=1/2 spins with four
% types of exchange coupling tensor: isotropic, two anisotropic, and
% antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-
% phonon relaxation in the generalised Lindblad form of Saito and
% Miyashita. The out-of-equilibrium curves are compared with the
% thermal equilibrium magnetisation. Reproduces Figure 4 of
%
%         https://arxiv.org/abs/2609.16352
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function dimer_exchange_types()

% Magnet must be 1 Tesla, the field is set by the sweep
sys.magnet=1.0;

% Parallel pool size
sys.parallel={'processes',4};

% Two electron spins with g=2
sys.isotopes={'E','E'};
inter.zeeman.scalar={2.0 2.0};

% Temperature of the phonon bath
inter.temperature=0.2;

% Formalism and basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Exchange coupling tensors of the paper, cm^-1, in its H=-2*S1*J*S2 convention
tensors={[0.2 0 0; 0 0.2 0; 0 0 0.2], ...
         [0.2 0 0; 0 0   0; 0 0 0  ], ...
         [0   0 0; 0 0   0; 0 0 0.2], ...
         [0 0.2 0.2; -0.2 0 0.2; -0.2 -0.2 0]};
labels={'isotropic','anisotropic, $J_{xx}$','anisotropic, $J_{zz}$','antisymmetric'};

% Spin-phonon coupling operator: unit elements between adjacent total S_z states
Sz=kron(full(stevens(2,1,0)),eye(2))+kron(eye(2),full(stevens(2,1,0))); msz=diag(Sz);
parameters.phonon_x=double(abs(msz-msz.')==1);

% Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-10 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-10*1e12*(1e-12)^2*0.1883651568463003^2;

% Observable: total magnetic moment along Z in Bohr magnetons
parameters.coil=-2.0*Sz;

% Sweep: 10 T/ms to 1 T in 10 ns stairs, output every 10 stairs
parameters.field_prof=@(t) 1e4*t;
parameters.timestep=1e-8; parameters.nsteps=1e4; parameters.nout=10;

% Single crystal in the frame of the exchange tensor
parameters.spins={'E'}; parameters.orientation=[0 0 0];
parameters.needs={'zeeman_op'};

% Loop over the exchange tensors
kfigure(); scale_figure([2.0 1.6]); answers=cell(1,4);
for n=1:4

    % Spinach coupling convention is S1*A*S2 with A in Hz
    inter.coupling.matrix=cell(2,2);
    inter.coupling.matrix{1,2}=-2*icm2hz(tensors{n});

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Run the simulation
    answers{n}=crystal(spin_system,@pulsed_field,parameters,'labframe');

    % Thermal equilibrium magnetisation at the same fields
    [I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0]);
    Z=hamiltonian(assume(spin_system,'labframe','zeeman')); H0=H0-Z; m_eq=zeros(size(answers{n}.field));
    for k=1:numel(m_eq)
        H=full(H0+answers{n}.field(k)*Z); H=(H+H')/2; [V,E]=eig(H,'vector');
        pops=exp(-spin_system.tols.hbar*(E-min(E))/(spin_system.tols.kbol*inter.temperature));
        m_eq(k)=real(trace(parameters.coil'*(V*diag(pops/sum(pops))*V')));
    end
    answers{n}.obs_eq=m_eq;

    % Plot the sweep and the equilibrium curves
    subplot(2,2,n); plot(answers{n}.field,answers{n}.obs); hold on;
    plot(answers{n}.field,m_eq,'--'); hold off; kgrid; xlim tight;
    kxlabel('Field, Tesla'); kylabel('Magnetisation, $\mu_B$'); ktitle(labels{n});
    klegend({'10 T/ms sweep','equilibrium'},'Location','northwest'); drawnow;

end

% Save the curves
save('dimer_exchange_types.mat','answers','labels');

end

