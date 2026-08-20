% Field dependence of geminate CIDNP from a radical pair in a viscous
% solvent, computed with Redfield theory and with the lifetime-shifted
% Nakajima-Zwanzig kernel. At a rotational correlation time of 1 ns and
% a singlet recombination rate of 1e8 Hz, the anisotropic hyperfine
% relaxation proceeds at a fair fraction of the recombination rate, and
% the pair drains before the bath decorrelates: the spectral densities
% seen by the surviving pair are lifetime-broadened, which changes the
% nuclear polarisation left in the diamagnetic product. The doubled-space
% bookkeeping follows the cidnp_geminate.m example; the low-field regime
% is motivated by the field-cycling CIDNP work of the Yurkovskaya and
% Ivanov school:
%
%          https://doi.org/10.1039/c3cp44098k
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function cidnp_nz_1()

% Field grid, recombination rate, and correlation time
fields=[0.002 0.005 0.01 0.02 0.035 0.05 0.075];
k_rec=3e8; tau_c=1e-9;

% Preallocate product polarisations
pol_rf=zeros(size(fields));
pol_nz=zeros(size(fields));

% Loop over the field grid
for n=1:numel(fields)

    % Compute both theories at this field
    pol_rf(n)=geminate_pol(fields(n),k_rec,tau_c,'redfield');
    pol_nz(n)=geminate_pol(fields(n),k_rec,tau_c,'naka-zwan');

end

% Report the two curves
disp('field (T), Redfield product polarisation, Nakajima-Zwanzig product polarisation:');
disp(num2str([fields' pol_rf' pol_nz'],'%14.6g'));

% Plot the two curves
kfigure(); plot(1e3*fields,[pol_rf' pol_nz'],'-o');
klegend({'Redfield','Nakajima-Zwanzig'},'Location','southeast');
kxlabel('magnetic field, mT'); kgrid;
kylabel('product nuclear polarisation');

end

% Nuclear polarisation left in the product by the geminate stage
function pol=geminate_pol(field,k_rec,tau_c,theory)

% Magnet and isotopes
sys.magnet=field;
sys.isotopes={'E','E','1H'};

% Silence the console
sys.output='hush';

% Electron g-factors
inter.zeeman.scalar={2.0023 2.0034 0};

% Anisotropic hyperfine tensor on the proton
inter.coupling.eigs=cell(3,3);
inter.coupling.euler=cell(3,3);
inter.coupling.eigs{2,3}=mt2hz([0.6 0.6 3.6]);
inter.coupling.euler{2,3}=[0 0 0];

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

% Singlet-born initial condition
rho=singlet(spin_system,1,2);

% Double up the problem, products carry no dynamics
H=[1*H 0*H; 0*H 0*H];
R=[1*R 0*R; 0*R 0*R];
rho=[1*rho; 0*rho];

% Reaction moves whatever leaves the reactants into the products
K=[1*K 0*K; -1*K 0*K];

% Assemble the Liouvillian
L=H+1i*R+1i*K;

% Evolve through the geminate stage
rho=evolution(spin_system,L,[],rho,200e-9,1,'final');

% Nuclear magnetisation in the product block
Nz=state(spin_system,'Lz','1H');
rho_prod=rho((numel(rho)/2+1):end);
pol=real(Nz'*rho_prod);

end

