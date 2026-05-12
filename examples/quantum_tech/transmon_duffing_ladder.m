% Duffing-model energy ladder of a weakly anharmonic transmon,
% showing how the transition frequencies separate as anharmo-
% nicity increases. Inspired by the transmon model of Koch et
% al., Phys. Rev. A 76, 042319 (2007).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_duffing_ladder()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T5'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Transmon number and anharmonicity operators
N=operator(spin_system,'N',1);
K=operator(spin_system,'CCAA',1);

% Model parameters
omega=2*pi*5.0e9;
anharm=2*pi*linspace(-400e6,-50e6,80);

% Preallocate transition frequency array
trans_frq=zeros(numel(anharm),4);

% Sweep the Duffing anharmonicity
for n=1:numel(anharm)

    % Build the Duffing Hamiltonian
    H=omega*N+(anharm(n)/2)*K;

    % Extract adjacent transition frequencies
    energies=sort(real(eig(full(H))));
    trans_frq(n,:)=diff(energies)/(2*pi);

end

% Plot the transition ladder
kfigure(); plot(-anharm/(2*pi*1e6),trans_frq/1e9,'LineWidth',1.5);
axis tight; kgrid; kxlabel('$-\alpha/2\pi$, MHz');
kylabel('transition frequency, GHz');
ktitle('Duffing transmon ladder');
klegend({'0-1','1-2','2-3','3-4'},'Location','SouthWest');

end

