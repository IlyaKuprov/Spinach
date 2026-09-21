% Pulsed-field magnetisation of a dimer of two S=1/2 spins with four
% types of exchange coupling tensor: isotropic, two anisotropic, and
% antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-
% phonon relaxation in the generalised Lindblad form of Saito and
% Miyashita. The out-of-equilibrium curves are compared with the
% thermal equilibrium magnetisation. Reproduces Fig 4 of
%
%                  https://arxiv.org/abs/2609.16352
%
% with the colours, line styles, and axis limits of that paper.
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

% Curve colours and legend labels of the paper
colours={[1 0 0],[0 0 1],[0 0.6 0],[1 0.65 0]}; 
labels=cell(1,8);

% Super-Ohmic bath, lambda^2*I0 of the paper 
% (lambda=10 cm^-1, I0=1e-10 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-10*1e12*(1e-12)^2*0.1883651568463003^2;

% Sweep: 10 T/ms to 1 T in 10 ns 
% stairs, output every 10 stairs
parameters.field_prof=@(t) 1e4*t;
parameters.timestep=1e-8; parameters.nsteps=1e4; parameters.nout=10;

% Single crystal in the frame of the exchange tensor
parameters.spins={'E'}; parameters.orientation=[0 0 0];
parameters.needs={'zeeman_op'};

% Loop over the exchange tensors
kfigure(); hold on; answers=cell(1,4);
for n=1:4

    % Spinach coupling convention is S1*A*S2 with A in Hz
    inter.coupling.matrix=cell(2,2);
    inter.coupling.matrix{1,2}=-2*icm2hz(tensors{n});

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Spin-phonon coupling operator: unit elements 
    % between adjacent total S_z states, projections 
    % rounded to exact integers
    Sz=full(operator(spin_system,'Lz','E')); 
    msz=round(2*diag(Sz))/2;
    parameters.phonon_x=double(abs(msz-msz.')==1);

    % Observable: total magnetic moment 
    % along Z in Bohr magnetons
    parameters.coil=-2.0*Sz;

    % Run the simulation
    answers{n}=crystal(spin_system,@pulsed_field,parameters,'labframe');

    % Thermal equilibrium magnetisation at the same fields
    [I,Q]=hamiltonian(assume(spin_system,'labframe')); 
    H0=I+orientation(Q,[0 0 0]);
    Z=hamiltonian(assume(spin_system,'labframe','zeeman')); 
    H0=H0-Z; m_eq=zeros(size(answers{n}.field));
    for k=1:numel(m_eq)
        H=full(H0+answers{n}.field(k)*Z); 
        rho=equilibrium(spin_system,(H+H')/2);
        m_eq(k)=real(hdot(parameters.coil,rho));
    end
    answers{n}.obs_eq=m_eq;

    % Equilibrium as a solid line and the sweep as a
    % dashed line in the colour of the tensor
    plot(answers{n}.field,m_eq,'-','Color',colours{n},'LineWidth',1.5);
    plot(answers{n}.field,answers{n}.obs,'--',...
         'Color',colours{n},'LineWidth',1.5);
    labels{2*n-1}=['$\mathbf{J}_' num2str(n) '^{\rm dimer}$ Equilibrium'];
    labels{2*n}=['$\mathbf{J}_' num2str(n) '^{\rm dimer}$ QME']; drawnow;

end

% Axis limits and legend of the paper
hold off; kgrid; xlim([0 1]); xticks(0:0.25:1);
ylim([0 2]); yticks(0:0.5:2);
kxlabel('$B$ (T)'); kylabel('Magnetisation ($\mu_B$)');
klegend(labels,'Location','east');

end

