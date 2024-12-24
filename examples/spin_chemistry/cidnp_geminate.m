% A basic example of the geminate CIDNP effect simulation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function cidnp_geminate()

% System specification
sys.magnet=14.1;
sys.isotopes ={'E','E','1H'};
inter.zeeman.scalar={2.0023 2.0024 1.0};
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{2,3}=1e7;
inter.chem.rp_theory='haberkorn';
inter.chem.rp_electrons=[1 2];
inter.chem.rp_rates=[1e7 0];
                     
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the Hamiltonian
H=hamiltonian(assume(spin_system,'esr'));

% Get the kinetics superoperator
K=kinetics(spin_system);

% Get the initial state
rho=singlet(spin_system,1,2);

% Double up the problem (no dynamics assumed in the product subspace)
H=[1*H 0*H; 0*H 0*H];
rho=[1*rho; 0*rho];

% Set up a reaction ("whatever is leaving reactants must appear in products")
K=[1*K 0*K; -1*K  0*K];

% Assemble the Liouvillian
L=H+1i*K;

% Evolve for a microsecond
rho=evolution(spin_system,L,[],rho,1e-6,1,'final');

% Check the nuclear magnetisation
Nz=state(spin_system,'Lz','1H');
rho_reac=rho(1:(numel(rho)/2));
rho_prod=rho((numel(rho)/2+1):end);
disp(['Nuclear magnetisation in reactants: ' num2str(real(Nz'*rho_reac))]);
disp(['Nuclear magnetisation in products:  ' num2str(real(Nz'*rho_prod))]);

end

