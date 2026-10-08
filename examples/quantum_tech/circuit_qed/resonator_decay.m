% Open-system dynamics of a leaky microwave resonator at finite
% temperature. A Fock state decays as a downward cascade through
% the level ladder, and a coherent state decays with its Poisson
% population structure largely preserved; both settle into the
% thermal state of the mode. The mean photon number follows the
% analytical amplitude damping solution in both cases, up to the
% distortion of the weak thermal channel by the Fock space trun-
% cation. Model and parameters from the resonator decay example
% of the paraqeet package.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function resonator_decay()

% Magnet field
sys.magnet=0;

% Microwave resonator with five Fock levels
sys.isotopes={'C5'};

% Resonator frequency
inter.modes.frqs={6.02e9};

% Energy relaxation lifetime
inter.modes.lifetimes={10e-9};

% Temperature of the environment
inter.temperature=0.050;

% Bose-Einstein occupation at the environment temperature
n_eq=1/(exp(6.62607015e-34*inter.modes.frqs{1}/...
            (1.380649e-23*inter.temperature))-1);

% Total coherence time from a five microsecond pure dephasing time
inter.modes.t2_times={1/(1/5e-6+(1+2*n_eq)/(2*inter.modes.lifetimes{1}))};

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Dissipative evolution generator
H=hamiltonian(assume(spin_system,'labframe'));
R=relaxation(spin_system);
G=-1i*H+R;

% Fock state and coherent state as initial conditions
rho_fock=state(spin_system,'BL5',1);
rho_coh=coherent(spin_system,1,1.5);

% Fock level population detection states
coils=[state(spin_system,'BL1',1) state(spin_system,'BL2',1) ...
       state(spin_system,'BL3',1) state(spin_system,'BL4',1) ...
       state(spin_system,'BL5',1)];

% Time grid of the source calculation
nsteps=100; dt=1e-9; time_axis=dt*(0:nsteps);

% Propagator over one time step
P=expm(full(G)*dt);

% Preallocate population trajectories
pops_fock=zeros(5,nsteps+1); pops_coh=zeros(5,nsteps+1);

% Propagate and record Fock level populations
for n=1:(nsteps+1)
    pops_fock(:,n)=real(coils'*rho_fock);
    pops_coh(:,n)=real(coils'*rho_coh);
    rho_fock=P*rho_fock; rho_coh=P*rho_coh;
end

% Mean photon numbers from the populations
n_fock=(0:4)*pops_fock; n_coh=(0:4)*pops_coh;

% Validate against the analytical solution, Fock truncation limits accuracy
kappa=1/inter.modes.lifetimes{1};
n_fock_ref=(n_fock(1)-n_eq)*exp(-kappa*time_axis)+n_eq;
n_coh_ref=(n_coh(1)-n_eq)*exp(-kappa*time_axis)+n_eq;
if max(abs(n_fock-n_fock_ref))>5e-3
    error('Fock state photon number deviates from the analytical solution.');
end
if max(abs(n_coh-n_coh_ref))>5e-3
    error('coherent state photon number deviates from the analytical solution.');
end

% Validate the downward cascade order of the population maxima
[~,idx_three]=max(pops_fock(4,:)); [~,idx_two]=max(pops_fock(3,:));
[~,idx_one]=max(pops_fock(2,:));
if ~((idx_three<idx_two)&&(idx_two<idx_one))
    error('Fock cascade maxima are out of order.');
end

% Validate the approach to the thermal state
if abs(pops_fock(1,end)-1/(1+n_eq))>2e-3
    error('final state is not thermal.');
end

% Plot the population dynamics for both initial states
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(1e9*time_axis,pops_fock','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns'); kylabel('Fock level populations');
ktitle('Fock state $|4\rangle$ decay');
klegend({'$|0\rangle$','$|1\rangle$','$|2\rangle$',...
         '$|3\rangle$','$|4\rangle$'},'Location','Best');
subplot(1,2,2); plot(1e9*time_axis,pops_coh','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns'); kylabel('Fock level populations');
ktitle('coherent state $|\alpha=1.5\rangle$ decay');
klegend({'$|0\rangle$','$|1\rangle$','$|2\rangle$',...
         '$|3\rangle$','$|4\rangle$'},'Location','Best');

end

