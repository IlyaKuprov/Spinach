% Demonstrate the use of the auxiliary matrix algorithm in generating a
% gradient sandwich multiple-quantum filter. For further details see:
%
%            http://dx.doi.org/10.1016/j.jmr.2014.01.011
%
% ledwards@cbs.mpg.de

function pfg_test_2()

% Magnet and isotopes
sys.magnet=5.9;
sys.isotopes={'1H','1H','1H'};

% Random chemical shifts and couplings
inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)};
inter.coupling.scalar=cell(3);
inter.coupling.scalar{1,2}=20*rand(1);
inter.coupling.scalar{1,3}=20*rand(1);
inter.coupling.scalar{2,3}=20*rand(1);

% Set the basis 
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');

% Build the Hamiltonian
H=hamiltonian(spin_system);

% Propagator for a pi/2 pulse
Hp=operator(spin_system,'L+','1H');
P=propagator(spin_system,0.5*(Hp+Hp'),pi/2);

% Build initial state vector
rho=rand(size(H,1),1);
rho(1)=1; rho=rho./norm(rho);
 
% Determine projection quantum numbers of the basis
[~,M]=lin2lm(spin_system.bas.basis);

% Determine the coherence order of each state
coherence_orders=sum(M,2);

% Find out which coherence orders are present
unique_coherence_orders=unique(coherence_orders).';

% Choose initial state weightings
weighting=[0,0,0,0,0,1,0];

% Weight coherence orders by the number of states
weight_idx=1;
for coherence_order=unique_coherence_orders
    subspace_mask=(coherence_order==coherence_orders);
    rho(subspace_mask)=rho(subspace_mask)./norm(rho(subspace_mask))*weighting(weight_idx);
    weight_idx=weight_idx+1;
end

% Build coherence diagram
gradient_strengths=20*[1 1];  % G/cm
sample_length=2.5;            % cm
gradient_shape_factors=[1 1]; % Rectangular shape function
nsteps=1000; durs=[1 1]*1e-3;
time_axis=linspace(0,sum(durs),nsteps);

% Preallocate the trajectory
rho_stack=zeros([length(rho) nsteps],'like',1i);

% Compute the trajectory
parfor time_idx=1:nsteps
    durations_temp=durs;
    durations_temp(2)=time_axis(time_idx);  
    rho_stack(:,time_idx)=grad_sandw(spin_system,H,rho,P,gradient_strengths,sample_length,durations_temp,gradient_shape_factors);
end

% Analyze the trajectory
figure();
trajan(spin_system,rho_stack,'coherence_order');
set(gca,'yscale','linear');

end

