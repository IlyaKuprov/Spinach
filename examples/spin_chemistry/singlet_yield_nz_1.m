% Magnetic field effect on a triplet-born benzophenone ketyl / thiyl
% radical pair in a viscous ionic liquid, computed with Redfield theory
% and with the lifetime-shifted Nakajima-Zwanzig kernel. The cage life-
% time of tens of nanoseconds and the rotational correlation time of a
% few nanoseconds put k*tau_c near 0.1, where the recombination drain
% visibly competes with rotational decorrelation of the anisotropic
% hyperfine couplings. Parameters are representative of the TMPA-TFSA
% measurements of the Wakasa group:
%
%          https://doi.org/10.1021/jp074331a
%          https://doi.org/10.1039/c2cp23747d
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function singlet_yield_nz_1()

% Field grid, cage drain, and correlation time
fields=linspace(0.5e-3,50e-3,15);
k_cage=3e7; tau_c=3.3e-9;

% Preallocate singlet yields
yield_rf=zeros(size(fields));
yield_nz=zeros(size(fields));

% Loop over the field grid
for n=1:numel(fields)

    % Compute both theories at this field
    yield_rf(n)=cage_pair(fields(n),k_cage,tau_c,'redfield');
    yield_nz(n)=cage_pair(fields(n),k_cage,tau_c,'naka-zwan');

end

% Report the two curves
disp('field (T), Redfield yield, Nakajima-Zwanzig yield:');
disp(num2str([fields' yield_rf' yield_nz'],'%14.6g'));

% Drain rate ladder at the rising edge of the MARY curve
drains=[1e6 3e6 1e7 3e7 1e8];
edge_rf=zeros(size(drains)); edge_nz=zeros(size(drains));
for n=1:numel(drains)
    edge_rf(n)=cage_pair(3e-3,drains(n),tau_c,'redfield');
    edge_nz(n)=cage_pair(3e-3,drains(n),tau_c,'naka-zwan');
end

% Report the drain ladder
disp('drain rate (Hz), Redfield yield, Nakajima-Zwanzig yield:');
disp(num2str([drains' edge_rf' edge_nz'],'%14.6g'));

% Plot the MARY curves and the drain ladder
kfigure(); scale_figure([1.75 0.75]);
subplot(1,2,1); plot(1e3*fields,[yield_rf' yield_nz']);
klegend({'Redfield','Nakajima-Zwanzig'},'Location','northeast');
kxlabel('magnetic field, mT'); kgrid;
kylabel('singlet recombination yield');
subplot(1,2,2); semilogx(tau_c*drains,edge_nz-edge_rf,'-o');
kxlabel('$k\tau_{c}$'); kgrid;
kylabel('NZ minus Redfield yield');

end

% Singlet recombination yield of the caged pair
function yield=cage_pair(field,k_cage,tau_c,theory)

% Magnet and isotopes
sys.magnet=field;
sys.isotopes={'E','E','1H','1H'};

% Silence the console
sys.output='hush';

% Electron g-factors, ketyl and thiyl
inter.zeeman.scalar={2.0032 2.0080 0 0};

% Anisotropic ketyl proton hyperfine tensors
inter.coupling.eigs=cell(4,4);
inter.coupling.euler=cell(4,4);
inter.coupling.eigs{1,3}=mt2hz([0.14 0.14 0.59]);
inter.coupling.eigs{1,4}=mt2hz([0.21 0.21 0.66]);
inter.coupling.euler{1,3}=[0 0 0];
inter.coupling.euler{1,4}=[pi/7 pi/5 pi/3];

% Exponential cage recombination kinetics
inter.chem.rp_theory='exponential';
inter.chem.rp_electrons=[1 2];
inter.chem.rp_rates=[k_cage/2 k_cage/2];

% Relaxation theory settings
inter.relaxation={theory};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='keep';
inter.tau_c={tau_c};

% Lifetime shift from the declared kinetics
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

% Singlet projector and triplet-born initial state
S=singlet(spin_system,1,2);
EE=state(spin_system,{'E','E'},{1,2});
rho0=(EE-S)/3;

% Resolvent solve for the reaction integral
L=H+1i*R+1i*K;
x=bicg(L,rho0,1e-8,numel(rho0));

% Convention-free fractional singlet yield
yield=real((k_cage/2)*imag(S'*x)/(k_cage*imag(EE'*x)));

end

