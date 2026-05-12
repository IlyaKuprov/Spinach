% Collective normal-mode splitting in the Tavis-Cummings model
% for one to four identical electron spins coupled to a common
% microwave cavity mode. The bright-state splitting follows the
% square-root scaling of Tavis and Cummings, Phys. Rev. 170,
% 379 (1968).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function tavis_cummings_splitting()

% Coupling strength
coupling=2*pi*3e6;

% Preallocate the splitting array
spin_counts=1:4;
splitting=zeros(size(spin_counts));

% Loop over ensemble sizes
for n=spin_counts

    % Magnet field
    sys.magnet=0;

    % Particle specification
    sys.isotopes=[repmat({'E'},1,n) {'C3'}];

    % Formalism and basis
    bas.formalism='zeeman-hilb';
    bas.approximation='none';

    % Spinach housekeeping
    spin_system=create(sys,[]);
    spin_system=basis(spin_system,bas);

    % Build the Tavis-Cummings Hamiltonian
    H=sparse(prod(spin_system.comp.mults),prod(spin_system.comp.mults));
    for k=1:n
        H=H+coupling*(operator(spin_system,{'L+','A'},{k,n+1})+...
                      operator(spin_system,{'L-','C'},{k,n+1}));
    end

    % Clean up numerical asymmetry
    H=(H+H')/2;

    % Extract the bright-mode splitting
    energies=sort(real(eig(full(H))));
    positive=min(energies(energies>1e-6));
    negative=max(energies(energies<-1e-6));
    splitting(n)=positive-negative;

end

% Plot numerical and analytical splittings
kfigure(); plot(spin_counts,splitting/(2*pi*1e6),'ko','MarkerFaceColor','k');
hold on; plot(spin_counts,2*sqrt(spin_counts)*coupling/(2*pi*1e6),'r-','LineWidth',1.5);
hold off; axis tight; kgrid; kxlabel('number of spins');
kylabel('normal-mode splitting, MHz');
ktitle('Tavis-Cummings square-root scaling');
klegend({'Spinach','2g\surd N'},'Location','NorthWest');

end

