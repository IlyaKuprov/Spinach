% Energy level diagram transition from the Zeeman limit to
% the exchange coupling limit in a two-electron system.
%
% i.kuprov@soton.ac.uk

function energy_levels()

% 600 MHz magnet
sys.magnet=14.1;

% Two electrons
sys.isotopes={'E','E'};

% Exaggerate g-factor difference
inter.zeeman.scalar={1.9 2.1};

% Hilbert space calculation
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relevant operators
Hz=hamiltonian(assume(spin_system,'labframe'));
Lj=operator(spin_system,{'Lx','Lx'},{1,2})+...
   operator(spin_system,{'Ly','Ly'},{1,2})+...
   operator(spin_system,{'Lz','Lz'},{1,2});

% omega_j scan from 0 to omega_e
omega_e=sys.magnet*spin('E');
omega_j=linspace(-3*omega_e,3*omega_e,1000);

% Get the energy levels
for n=1:numel(omega_j)
    
    % Diagonalise the Hamiltonian
    [~,D]=eig(Hz-omega_j(n)*Lj,'vector');
    
    % Sort and record energies
    energies(:,n)=sort(D); %#ok<AGROW>
    
end

% Plot energies
figure(); plot(omega_j/omega_e,energies'/omega_e,'k-'); 
kxlabel('$\omega_{J}/\omega_{E}$');
kylabel('$\omega/\omega_{E}$'); kgrid;

end

