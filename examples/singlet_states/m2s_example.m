% An example of the M2S sequence for a two-spin system.
%
% Calculation time: seconds
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function m2s_example()

% Spin system and interactions
sys.magnet=9.4;
sys.isotopes={'13C','13C'};
inter.zeeman.scalar={0.03,-0.03};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=55;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Pulse operators
Hy=(operator(spin_system,'L+','13C')-...
    operator(spin_system,'L-','13C'))/2i;
Hx=(operator(spin_system,'L+','13C')+...
    operator(spin_system,'L-','13C'))/2;

% Start with longitudinal magnetisation
rho0=state(spin_system,'Lz','all');

% Detect singlet state
coil=singlet(spin_system,1,2);

% Call the M2S sequence
rho=m2s(spin_system,H,Hx,Hy,rho0,55,6.0);

% Display the singlet population
disp(['Singlet population: ' num2str(coil'*rho)]);

end

