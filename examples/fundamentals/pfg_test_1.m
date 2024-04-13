% A test of the explicit gradient pulse function that uses the auxiliary 
% matrix formalism to compute sample volume integral. For details, see:
%
%               http://dx.doi.org/10.1016/j.jmr.2014.01.011
%
% luke.edwards@chem.ox.ac.uk

function pfg_test_1()

% Magnet and isotopes
sys.magnet=5.9;
sys.isotopes={'1H','1H','1H'};

% Random chemical shifts and couplings
inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)};
inter.coupling.scalar=cell(3);
inter.coupling.scalar{1,2}=20*rand(1);
inter.coupling.scalar{1,3}=20*rand(1);
inter.coupling.scalar{2,3}=20*rand(1);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');

% Spin Hamiltonian
L=hamiltonian(spin_system);

% Build initial state vector
rho=rand(size(L,1),1);
rho(1)=1; rho=rho./norm(rho);
 
% Determine projection quantum numbers of the basis
[~,M]=lin2lm(spin_system.bas.basis);

% Determine the coherence order of each state
coherence_orders=sum(M,2);

% Find out which coherence orders are present
unique_coherence_orders=unique(coherence_orders).';

% Weight coherence orders by the number of states
weighting=linspace(0.1,0.9,numel(unique_coherence_orders)); weight_idx=1;
for coherence_order=unique_coherence_orders
    subspace_mask=(coherence_order==coherence_orders);
    rho(subspace_mask)=rho(subspace_mask)./norm(rho(subspace_mask))*weighting(weight_idx);
    weight_idx=weight_idx+1;
end

% Evolve under the homospoil pulse
nsteps=100;
delta_step=2e-7; % s
gradient_strength=20; % G/cm
sample_length=1.5; % cm
gradient_shape_factor=1; % Rectangular shape function
time_axis=0:delta_step:(nsteps-1)*delta_step;

% Preallocate the trajectory
rho_stack=zeros([length(rho) nsteps],'like',1i);

% Compute the trajectory
parfor t_idx=1:nsteps
    duration=time_axis(t_idx);
    rho_stack(:,t_idx)=grad_pulse(spin_system,L,rho,gradient_strength,sample_length,duration,gradient_shape_factor);
end

% Analyze the trajectory
figure();
trajan(spin_system,rho_stack,'coherence_order');
set(gca,'yscale','linear');

end

