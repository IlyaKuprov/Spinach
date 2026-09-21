% Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a
% linear trimer of three S=5/2 manganese ions with isotropic exchange
% between neighbours and an axial plus rhombic zero-field splitting on
% every ion, from zero to 10 Tesla in the full 216-state Hilbert space.
% The lowest levels are the ones that the 16-state and 26-state effec-
% tive bases of the paper are built to reproduce. Reproduces Fig 5 of
%
%                  https://arxiv.org/abs/2609.16352
%
% with the axis limits, colours, and the slight vertical offsets that
% the paper applies to show the overlapping lines; the 16-state and
% 26-state effective bases of the paper (common eigenstates of the
% isotropic exchange and the total S_z) come from mn3_trimer_basis.m.
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

% Isotropic exchange, J=-2.42 cm^-1 in 
% the H=-2*J*S1*S2 convention of the paper
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

% Field-free Hamiltonian and the Zeeman operator per Tesla
[I,Q]=hamiltonian(assume(spin_system,'labframe')); 
Z=hamiltonian(assume(spin_system,'labframe','zeeman')); 
H0=I+orientation(Q,[0 0 0]); H0=H0-Z;

% Full space and the two effective bases of the paper
bases={speye(216),mn3_trimer_basis(spin_system,16),...
       mn3_trimer_basis(spin_system,26)};

% Energy levels on a field grid, cm^-1, in each basis
fields=linspace(0,10,201); levels=cell(1,3);
for b=1:3
    levels{b}=zeros(size(bases{b},2),numel(fields));
    for k=1:numel(fields)
        H=full(bases{b}'*(H0+fields(k)*Z)*bases{b});
        levels{b}(:,k)=hz2icm(sort(eig((H+H')/2))/(2*pi));
    end
end

% Levels relative to the field-free ground state, offset
% by 0.15 cm^-1 per basis to show the overlapping lines
kfigure(); hold on; handles=zeros(1,3);
colours={[0.8 0.8 0.8],[0 0 0],[1 0 0]};
for b=1:3
    h=plot(fields,levels{b}-levels{1}(1,1)+0.15*(b-1),'-','Color',colours{b});
    handles(b)=h(1);
end
hold off; kgrid; xlim([0 10]); xticks(0:2:10);
ylim([-10 10]); yticks(-10:5:10);
kxlabel('$B$ (T)'); kylabel('$E$ (cm$^{-1}$)');
klegend(handles,{'All states','16 states','26 states'},...
        'Location','southwest');

end

