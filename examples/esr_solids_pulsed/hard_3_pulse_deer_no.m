% Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals
% are positioned at a distance of 25 Angstroms.
%
% The numerical calculation is done by brute-force time propaga-
% tion and numerical powder averaging, including g-factor orien-
% tation effects on the dipolar coupling. Nitroxide g-tensor is
% from http://dx.doi.org/10.1063/1.1697233
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function hard_3_pulse_deer_no()

% Spin system properties
sys.magnet=0.33;
sys.isotopes={'E','E'};
inter.zeeman.eigs{1}=[2.0089 2.0061 2.0027];
inter.zeeman.eigs{2}=[2.0089 2.0061 2.0027];
inter.zeeman.euler{1}=[1.0 2.0 3.0];
inter.zeeman.euler{2}=[3.0 1.0 2.0];
inter.coordinates={[ 0.00 0.00 0.00]
                   [25.00 0.00 0.00]};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';
               
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil_prob=state(spin_system,{'L-'},{1});
parameters.stepsize=1e-8;
parameters.nsteps=100;
parameters.spins={'E'};
parameters.ex_prob=(operator(spin_system,{'L+'},{1})+...
                    operator(spin_system,{'L-'},{1}))/2;
parameters.ex_pump=(operator(spin_system,{'L+'},{2})+...
                    operator(spin_system,{'L-'},{2}))/2;
parameters.output='brief';
parameters.grid='rep_2ang_3200pts_sph';

% Pulse sequence
deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer');

% Build the time axis
time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1);

% Plotting
figure(); plot(1e6*time_axis,imag(deer.deer_trace)); 
xlabel('time, microseconds'); axis tight; kgrid; 

end

