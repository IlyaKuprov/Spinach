% Redfield superoperator for the scalar relaxation of the first
% kind in a two-proton system with a noisy J-coupling. This si-
% tuation occurs in aziridines, where the slow nitrogen inver-
% sion jitters scalar couplings on a millisecond time scale.
% Set to demonstrate the effect described in:
%
%           http://dx.doi.org/10.1002/ange.201410271
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% tim.claridge@chem.ox.ac.uk
% barbara.odell@chem.ox.ac.uk

function scalar_relaxation_1()

% System specification
sys.magnet=11.75;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0.0 2.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation superoperator
inter.relaxation={'SRFK'};
inter.rlx_keep='kite';
inter.equilibrium='zero';
inter.srfk_tau_c={[1.0 1e-3]};
inter.srfk_mdepth=cell(2);
inter.srfk_mdepth{1,2}=15.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Show a spy plot of R
figure(); spy(relaxation(spin_system));

end

