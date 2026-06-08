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

    % Locate the one-excitation manifold
    one_quant=speye(size(H,1));
    subspace=zeros(n+1,1);
    for k=1:n
        spin_state=repmat({'ZL1'},1,n);
        spin_state{k}='ZL2';
        projector=state(spin_system,[spin_state {'BL1'}],num2cell(1:(n+1)));
        subspace(k)=find(diag(projector)>0.5);
    end
    projector=state(spin_system,[repmat({'ZL1'},1,n) {'BL2'}],num2cell(1:(n+1)));
    subspace(n+1)=find(diag(projector)>0.5);
    one_quant=one_quant(:,subspace);

    % Extract the bright-mode splitting
    energies=sort(real(eig(full(one_quant'*H*one_quant))));
    splitting(n)=max(energies)-min(energies);

end

% Compute the analytical splitting
analytical=2*sqrt(spin_counts)*coupling;

% Validate the square-root scaling
if max(abs(splitting-analytical)./analytical)>1e-10
    error('Tavis-Cummings splitting does not follow square-root scaling.');
end

% Plot numerical and analytical splittings
kfigure(); plot(spin_counts,splitting/(2*pi*1e6),'ko','MarkerFaceColor','k');
hold on; plot(spin_counts,analytical/(2*pi*1e6),'r-','LineWidth',1.5);
hold off; axis tight; kgrid; kxlabel('number of spins');
kylabel('normal-mode splitting, MHz');
ktitle('Tavis-Cummings square-root scaling');
klegend({'Spinach','2g\surd N'},'Location','Best');

end
