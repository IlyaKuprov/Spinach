% Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a
% linear trimer of three S=5/2 manganese ions with isotropic exchange
% between neighbours and an axial plus rhombic zero-field splitting on
% every ion, from zero to 10 Tesla in the full 216-state Hilbert space.
% The lowest levels are the ones that the 16-state and 26-state effec-
% tive bases of the paper are built to reproduce. Reproduces Fig 5 of
%
%                  https://arxiv.org/abs/2609.16352
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function mn3_trimer_levels()

% Magnet must be 1 Tesla, the field is set below
sys.magnet=1.0;

% Parallel pool size
sys.parallel={'processes',4};

% Three S=5/2 spins with g=2
sys.isotopes={'E6','E6','E6'};
inter.zeeman.scalar={2.0 2.0 2.0};

% Isotropic exchange, J=-2.42 cm^-1 in the H=-2*J*S1*S2 convention of the paper
inter.coupling.matrix=cell(3,3);
inter.coupling.matrix{1,2}=-2*icm2hz(-2.42)*eye(3);
inter.coupling.matrix{2,3}=-2*icm2hz(-2.42)*eye(3);

% Zero-field splitting, D=0.167 cm^-1 and E=0.040 cm^-1 on every ion
for n=1:3
    inter.coupling.matrix{n,n}=zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0);
end

% Formalism and basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Field-free Hamiltonian and the Zeeman operator per Tesla
[I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0]);
Z=hamiltonian(assume(spin_system,'labframe','zeeman')); H0=H0-Z;

% Energy levels on a field grid, cm^-1
fields=linspace(0,10,201); levels=zeros(size(H0,1),numel(fields));
for k=1:numel(fields)
    H=full(H0+fields(k)*Z); levels(:,k)=hz2icm(sort(eig((H+H')/2))/(2*pi));
end

% Plot the lowest thirty levels relative to the field-free ground state
kfigure(); plot(fields,levels(1:30,:)-levels(1,1)); kgrid; xlim tight; ylim padded;
kxlabel('Field, Tesla'); kylabel('Energy, cm$^{-1}$');

end

