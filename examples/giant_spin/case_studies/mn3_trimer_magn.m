% Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a
% linear trimer of three S=5/2 manganese ions with isotropic exchange
% between neighbours and an axial plus rhombic zero-field splitting on
% every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon
% relaxation in the generalised Lindblad form of Saito and Miyashita.
% The full 216-state Hilbert space is used; the paper solves the same
% problem in 16-state and 26-state effective bases. The thermal equili-
% brium magnetisation is plotted for comparison. Reproduces Fig 6 of
%
%                 https://arxiv.org/abs/2609.16352
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

% Total S_z and the spin-phonon coupling 
% operator between adjacent total S_z states, 
% projections rounded to exact half-integers
Sz=full(operator(spin_system,'Lz','E6')); 
msz=round(2*diag(Sz))/2;
parameters.phonon_x=double(abs(msz-msz.')==1);

% Super-Ohmic bath, lambda^2*I0 of the paper 
% (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2;

% Observable: total magnetic 
% moment along Z in Bohr magnetons
parameters.coil=-2.0*Sz;

% Sweep: 50 T/ms to 10 T in 20 ns 
% stairs, output every 100 stairs
parameters.field_prof=@(t)5e4*t;
parameters.timestep=2e-8; 
parameters.nsteps=1e4; 
parameters.nout=100;

% Single crystal in the frame of
% the zero-field splitting tensors
parameters.spins={'E6'}; 
parameters.orientation=[0 0 0];
parameters.needs={'zeeman_op'};

% Run the simulation
answer=crystal(spin_system,@pulsed_field,parameters,'labframe');

% Thermal equilibrium magnetisation at the same fields
[I,Q]=hamiltonian(assume(spin_system,'labframe'));
Z=hamiltonian(assume(spin_system,'labframe','zeeman')); 
H0=I+orientation(Q,[0 0 0]); H0=H0-Z; 
m_eq=zeros(size(answer.field));
for k=1:numel(m_eq)
    H=full(H0+answer.field(k)*Z); 
    rho=equilibrium(spin_system,(H+H')/2);
    m_eq(k)=real(hdot(parameters.coil,rho));
end

% Plot the sweep and the equilibrium curves
kfigure(); plot(answer.field,answer.obs); 
hold on; plot(answer.field,m_eq,'--'); hold off;
kgrid; xlim tight; kxlabel('Field, Tesla'); 
kylabel('Magnetisation, $\mu_B$');
klegend({'50 T/ms sweep, 216 states','equilibrium'},...
        'Location','northwest');

end

