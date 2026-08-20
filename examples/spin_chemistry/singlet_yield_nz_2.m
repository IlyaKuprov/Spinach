% Field dependence of the decay rate of a micelle-confined triplet-born
% benzophenone ketyl / alkyl radical pair, computed with Redfield theory
% and with the lifetime-shifted Nakajima-Zwanzig kernel. The high-field
% decay of micellar pairs is relaxation-controlled: the T+/- states drain
% into the reactive S/T0 subspace at rates set by spectral densities of
% the anisotropic hyperfine modulation. Contact recombination at 1e9 Hz
% against a supercage correlation time of 0.7 ns puts the scalar lifetime
% shift at k*tau_c near 0.35, and the two theories separate. The observed
% decay rate is the slowest eigenmode of the pair Liouvillian. Parameters
% are representative of the SDS supercage systems of Sakaguchi, Hayashi,
% and Nagakura:
%
%          https://doi.org/10.1246/bcsj.57.322
%          https://doi.org/10.1016/0009-2614(80)80325-2
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function singlet_yield_nz_2()

% Field grid, contact recombination rate, and correlation time
fields=[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.7 1.0 1.34];
k_rec=1e9; tau_c=0.7e-9;

% Preallocate observed decay rates
decay_rf=zeros(size(fields));
decay_nz=zeros(size(fields));

% Loop over the field grid
for n=1:numel(fields)

    % Compute both theories at this field
    decay_rf(n)=scrp_decay(fields(n),k_rec,tau_c,'redfield');
    decay_nz(n)=scrp_decay(fields(n),k_rec,tau_c,'naka-zwan');

end

% Report the two curves
disp('field (T), Redfield decay rate (Hz), Nakajima-Zwanzig decay rate (Hz):');
disp(num2str([fields' decay_rf' decay_nz'],'%14.6g'));

% Plot the two curves
kfigure(); semilogy(fields,[decay_rf' decay_nz'],'-o');
klegend({'Redfield','Nakajima-Zwanzig'},'Location','northeast');
kxlabel('magnetic field, Tesla'); kgrid;
kylabel('slowest pair decay rate, Hz');

end

% Slowest decay eigenmode of the confined pair
function decay=scrp_decay(field,k_rec,tau_c,theory)

% Magnet and isotopes
sys.magnet=field;
sys.isotopes={'E','E','1H','1H'};

% Silence the console
sys.output='hush';

% Electron g-factors, ketyl and alkyl
inter.zeeman.scalar={2.0031 2.0026 0 0};

% Anisotropic alkyl proton hyperfine tensors
inter.coupling.eigs=cell(4,4);
inter.coupling.euler=cell(4,4);
inter.coupling.eigs{2,3}=mt2hz([1.7 1.7 3.2]);
inter.coupling.eigs{2,4}=mt2hz([2.1 2.1 3.3]);
inter.coupling.euler{2,3}=[0 0 0];
inter.coupling.euler{2,4}=[pi/5 pi/3 pi/7];

% Haberkorn recombination, singlet channel only
inter.chem.rp_theory='haberkorn';
inter.chem.rp_electrons=[1 2];
inter.chem.rp_rates=[k_rec 0];

% Relaxation theory settings
inter.relaxation={theory};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='keep';
inter.tau_c={tau_c};

% Scalar lifetime shift from the declared kinetics
if strcmp(theory,'naka-zwan')
    inter.nz_shift='chem';
    inter.nz_onshell=false;
end

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian, relaxation, and kinetics superoperators
H=hamiltonian(assume(spin_system,'labframe'));
R=relaxation(spin_system);
K=kinetics(spin_system);

% Decay rates of the Liouvillian eigenmodes
L=full(H+1i*R+1i*K);
mode_rates=-imag(eig(L));

% Slowest decaying mode above the numerical floor
decay=min(mode_rates(mode_rates>1));

end

