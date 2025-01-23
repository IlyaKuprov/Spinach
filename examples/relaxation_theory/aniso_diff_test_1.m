% Relaxation superoperator calculation for a dipole-coupled two-spin
% system with an anisotropic rotational diffusion tensor.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function aniso_diff_test_1()

% Magnet field (Tesla)
sys.magnet=2*pi*950.33e6/spin('1H');

% Isotopes
sys.isotopes={'1H','13C'};

% Chemical shift tensors (ppm)
inter.zeeman.eigs={[0.0  0.0  0.0];
                   [0.0  0.0  0.0]};
inter.zeeman.euler={[0.0  0.0  0.0];
                    [0.0  0.0  0.0]}; 

% Cartesian coordinates (Angstrom)
R=euler2dcm(1,2,3);
inter.coordinates={[1.13 0.00 0.00]*R;
                   [0.00 0.00 0.00]};
               
% Scalar couplings (Hz)
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=145.0;
                
% Difusion tensor eigenvalues
D=[2.16e8 2.35e8 7.45e8];

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1./(6*D)};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);
disp(full(R));

end

