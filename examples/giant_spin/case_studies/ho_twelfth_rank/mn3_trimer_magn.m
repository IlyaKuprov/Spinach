% Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a
% linear trimer of three S=5/2 manganese ions with isotropic exchange
% between neighbours and an axial plus rhombic zero-field splitting on
% every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon
% relaxation in the generalised Lindblad form of Saito and Miyashita.
% As in the paper, the dynamics runs in the 16-state and 26-state
% effective bases of mn3_trimer_basis.m, with the Hamiltonian, the
% Zeeman operator, and the observable projected onto each basis;
% the thermal equilibrium magnetisation of the full 216-state space
% is plotted for comparison. Reproduces Fig 6 of
%
%                 https://arxiv.org/abs/2609.16352
%
% with the colours and axis limits of that paper.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function mn3_trimer_magn()

% Magnet must be 1 Tesla, the field is set by the sweep
sys.magnet=1.0;

% Three S=5/2 spins with g=2, 
% phonon bath at 0.6 K
sys.isotopes={'E6','E6','E6'};
inter.zeeman.scalar={2.0 2.0 2.0};
inter.temperature=0.6;

% Isotropic exchange, J=-2.42 cm^-1 
% in the H=-2*J*S1*S2 convention of the paper
inter.coupling.matrix=cell(3,3);
inter.coupling.matrix{1,2}=-2*icm2hz(-2.42)*eye(3);
inter.coupling.matrix{2,3}=-2*icm2hz(-2.42)*eye(3);

% Zero-field splitting, D=0.167 cm^-1 
% and E=0.040 cm^-1 on every ion
for n=1:3
    inter.coupling.matrix{n,n}=zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0);
end

% Formalism and basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Total S_z operator
Sz=full(operator(spin_system,'Lz','E6'));

% Super-Ohmic bath, lambda^2*I0 of the paper 
% (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2;

% Observable: total magnetic 
% moment along Z in Bohr magnetons
coil=-2.0*Sz;

% Sweep: 50 T/ms to 10 T in 20 ns 
% stairs, output every 100 stairs
parameters.field_prof=@(t)5e4*t;
parameters.timestep=2e-8; 
parameters.nsteps=1e4; 
parameters.nout=100;

% Full Hamiltonian at 1 Tesla in the frame of the
% zero-field splitting tensors, Zeeman operator per Tesla
spin_system=assume(spin_system,'labframe');
[I,Q]=hamiltonian(spin_system); H_full=I+orientation(Q,[0 0 0]);
Z=hamiltonian(assume(spin_system,'labframe','zeeman'));

% Loop over the 16-state and 26-state
% effective bases of the paper
nstates=[16 26]; answers=cell(1,2);
for n=1:2

    % Basis states and their total S_z projections
    [P,msz]=mn3_trimer_basis(spin_system,nstates(n));

    % Hamiltonian, Zeeman operator, and
    % observable projected onto the basis
    H=P'*H_full*P; H=full((H+H')/2);
    parameters.hzeeman=full(P'*Z*P);
    parameters.hzeeman=(parameters.hzeeman+parameters.hzeeman')/2;
    parameters.coil=P'*coil*P;
    parameters.coil=(parameters.coil+parameters.coil')/2;

    % Spin-phonon coupling operator: unit elements
    % between basis states with adjacent total S_z
    parameters.phonon_x=double(abs(msz-msz.')==1);

    % Run the simulation
    answers{n}=pulsed_field(spin_system,parameters,H,[],[]);

end

% Thermal equilibrium magnetisation of
% the full space at the same fields
H0=H_full-Z; m_eq=zeros(size(answers{1}.field));
for k=1:numel(m_eq)
    H=full(H0+answers{1}.field(k)*Z);
    rho=equilibrium(spin_system,(H+H')/2);
    m_eq(k)=real(hdot(coil,rho));
end

% Equilibrium in blue, the 16-state sweep in
% red, and the 26-state sweep in green
kfigure(); plot(answers{1}.field,m_eq,'b-','LineWidth',1.5);
hold on; plot(answers{1}.field,answers{1}.obs,'r-','LineWidth',1.5);
plot(answers{2}.field,answers{2}.obs,'-','Color',[0 0.6 0],'LineWidth',1.5);
hold off; kgrid; xlim([0 10]); xticks(0:2.5:10); ylim([0 6]);
kxlabel('$B$ (T)'); kylabel('Magnetisation ($\mu_B$)');
klegend({'Equilibrium','QME - 16 states','QME - 26 states'},...
        'Location','southeast');

end

