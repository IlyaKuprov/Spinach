% Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin 
% complex with a strong exchange coupling between the electrons. A
% distribution in the exchange coupling is summed over.
%
% The calculation is done by brute-force time propagation and
% numerical powder averaging in Liouville space.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% sabine.richert@chem.ox.ac.uk
% christiane.timmel@chem.ox.ac.uk

function hard_3_pulse_deer_exchange()

% Generate the distribution
j_values=linspace(6e6,20e6,20);
weights=gaussfun(j_values-13.1e6,4.2622e6);
weights=weights/sum(weights); answer=0;

% Run the averaging
for n=1:numel(weights)
    
    % Hush up
    sys.output='hush';
    
    % Magnet field
    sys.magnet=1.2132;
    
    % Isotopes
    sys.isotopes={'E','E'};
    
    % Zeeman interactions
    inter.zeeman.eigs=cell(1,2);
    inter.zeeman.euler=cell(1,2);
    inter.zeeman.eigs{1}=[2.050 2.050 2.195];
    inter.zeeman.euler{1}=[0 0 0];
    inter.zeeman.eigs{2}=[2.050 2.050 2.195];
    inter.zeeman.euler{2}=[0 0 0];
    
    % Exchange coupling
    inter.coupling.scalar{1,2}=j_values(n);
    inter.coupling.scalar{2,2}=[];
    
    % Coordinates
    inter.coordinates=cell(2,1);
    inter.coordinates{1}=[0.00  0.00 0.00];
    inter.coordinates{2}=[24.50 0.00 0.00];
    
    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='none';
    
    % Disable trajectory level SSR algorithms
    sys.disable={'hygiene','trajlevel'};
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Sequence parameters
    parameters.rho0=state(spin_system,'Lz','E');
    parameters.coil_prob=state(spin_system,{'L-'},{1});
    parameters.stepsize=1e-8/4;
    parameters.nsteps=4*50;
    parameters.spins={'E'};
    parameters.ex_prob=operator(spin_system,{'Lx'},{1});
    parameters.ex_pump=operator(spin_system,{'Lx'},{2});
    parameters.output='brief';
    parameters.grid='rep_2ang_400pts_sph';
    
    % Pulse sequence
    deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer-zz');
    
    % Add to the total
    answer=answer+weights(n)*deer.deer_trace;

end

% Build time axis
time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1);

% Plotting
figure(); plot(1e6*time_axis,imag(answer)); 
axis tight; kgrid; kxlabel('time, microseconds');

end

