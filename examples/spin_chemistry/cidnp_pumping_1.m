% A simulation of the matrix in Equation 2 of IK's paper on chemically
% amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function cidnp_pumping_1()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','19F'};

% Chemical shifts
inter.zeeman.scalar={0.0, 0.0};

% Chemical shift anisotropies (DFT)
inter.zeeman.eigs{1}=[0 0 0];
inter.zeeman.euler{1}=[0 0 0];
inter.zeeman.eigs{2}=[-47 -16  63];
inter.zeeman.euler{2}=[0 0 0];

% Coordinates (DFT)
inter.coordinates={[0 0 0],[0 2.60 0]};

% J-coupling (expt)
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=50;

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.tau_c={110e-12};
inter.equilibrium='IME';
inter.temperature=298;

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get relaxation matrix
R=relaxation(spin_system);

% Build relevant states
U=unit_state(spin_system);        
Hz=state(spin_system,{'Lz'},{1}); 
Fz=state(spin_system,{'Lz'},{2}); 
HzFz=state(spin_system,{'Lz','Lz'},{1,2});
U=U/norm(U,2); Hz=Hz/norm(Hz,2);
Fz=Fz/norm(Fz,2); HzFz=HzFz/norm(HzFz,2);

% Add pumping terms
R=magpump(spin_system,R,Hz,1.3);
R=magpump(spin_system,R,Fz,34.0);

% Build active space projector
P=[U Hz Fz -HzFz];

% Project the relaxation matrix into the active space
disp('Kuprov''s matrix in Equation 2:'); disp(P'*R*P);

end

