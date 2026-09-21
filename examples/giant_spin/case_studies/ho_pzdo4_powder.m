% Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal-
% organic framework, a J=8 giant spin with a crystal field to twelfth
% spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon
% relaxation in the generalised Lindblad form of Saito and Miyashita.
% The magnetisation of every orientation of the two-angle Lebedev grid
% is plotted alongside the powder average and the thermal equilibrium
% magnetisation. Reproduces Figure 3 of
%
%         https://arxiv.org/abs/2609.16352
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function ho_pzdo4_powder()

% Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
[ks,qs,bkq]=ho_pzdo4_params();

% Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank
coeff=cell(1,12); euler=cell(1,12);
for k=1:12
    stev=zeros(2*k+1,1); sel=(ks==k); stev(qs(sel)+k+1)=bkq(sel);
    coeff{k}=icm2hz(stev2sph(k,stev)); euler{k}=[0 0 0];
end

% Magnet must be 1 Tesla, the field is set by the sweep
sys.magnet=1.0;

% J=8 giant spin, effective g-factor 1.24, phonon bath at 2 K
sys.isotopes={'E17'};
inter.zeeman.scalar={1.24};
inter.giant.coeff={coeff};
inter.giant.euler={euler};
inter.temperature=2.0;

% Formalism and basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Spin-phonon coupling operator: unit elements between adjacent m_J states
Jz=full(stevens(17,1,0)); mj=diag(Jz);
parameters.phonon_x=double(abs(mj-mj.')==1);

% Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2;

% Observable: magnetic moment along Z in Bohr magnetons
parameters.coil=-1.24*Jz;

% Sweep: 10 T/ms to 10 T in 10 ns stairs, output every 1000 stairs
parameters.field_prof=@(t) 1e4*t;
parameters.timestep=1e-8; parameters.nsteps=1e5; parameters.nout=1000;

% Two-angle Lebedev grid, outputs of every orientation returned separately
parameters.spins={'E17'}; parameters.grid='leb_2ang_rank_29';
parameters.needs={'zeeman_op'}; parameters.sum_up=false;

% Run the simulation
[answers,sph_grid]=powder(spin_system,@pulsed_field,parameters,'labframe');

% Powder average of the magnetisation
fields=answers{1}.field; m_avg=zeros(size(fields));
for n=1:numel(answers)
    m_avg=m_avg+sph_grid.weights(n)*answers{n}.obs;
end

% Thermal equilibrium magnetisation, powder averaged, at the same fields
[I,Q]=hamiltonian(assume(spin_system,'labframe')); [ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'));
m_eq=zeros(size(fields));
for n=1:numel(answers)
    angles=[sph_grid.alphas(n) sph_grid.betas(n) sph_grid.gammas(n)];
    H0=I+orientation(Q,angles); Z=ZI+orientation(ZQ,angles); H0=H0-Z;
    for k=1:numel(fields)
        H=full(H0+fields(k)*Z); H=(H+H')/2; [V,E]=eig(H,'vector');
        pops=exp(-spin_system.tols.hbar*(E-min(E))/(spin_system.tols.kbol*inter.temperature));
        m_eq(k)=m_eq(k)+sph_grid.weights(n)*real(trace(parameters.coil'*(V*diag(pops/sum(pops))*V')));
    end
end

% Plot the single orientation curves, the powder average, and the equilibrium
kfigure(); hold on;
for n=1:numel(answers)
    plot(fields,answers{n}.obs,'Color',[0.8 0.8 0.8]);
end
plot(fields,m_avg,'LineWidth',2); plot(fields,m_eq,'--','LineWidth',2); hold off;
kgrid; xlim tight; kxlabel('Field, Tesla'); kylabel('Magnetisation, $\mu_B$');
klegend({'single orientations','powder average','equilibrium'},'Location','northwest');

% Save the curves
save('ho_pzdo4_powder.mat','fields','m_avg','m_eq','answers','sph_grid');

end

