% An example of the S2M sequence for a two-spin system.
%
% Calculation time: seconds
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function s2m_example()

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
Cx=operator(spin_system,'Lx','13C');
Cy=operator(spin_system,'Ly','13C');

% Start with singlet state
rho0=singlet(spin_system,1,2);

% Detect longitudinal magnetisation
coil=state(spin_system,'Lz','all');

% Call the S2M sequence
rho=s2m(spin_system,H,Cx,Cy,rho0,55,6.0);

% Display the longitudinal magnetisation
disp(['Longitudinal magnetisation: ' num2str(coil'*rho)]);

end

