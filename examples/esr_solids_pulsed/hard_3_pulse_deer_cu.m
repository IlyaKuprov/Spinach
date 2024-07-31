% Three-pulse DEER on a Cu(II)-NO two electron system at X-band.
%
% The numerical calculation is done by brute-force time propaga-
% tion and numerical powder averaging in Liouville space, inclu-
% ding g-factor orientation effects on the dipolar coupling.
%
% The analytical calculation is done for isotropic parts of the
% electron g-factors.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function hard_3_pulse_deer_cu()

% Spin system parameters
sys.magnet=0.33;
sys.isotopes={'E','E'};
inter.zeeman.eigs={[2.056, 2.056, 2.205];
                   [2.009, 2.006, 2.003]};
inter.zeeman.euler={[0 0 0]; [0 0 0]};
inter.coordinates={[0 0 0]; [20 0 0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory level SSR algorithms
sys.disable={'trajlevel'};
               
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil_prob=state(spin_system,{'L-'},{1});
parameters.stepsize=1e-8;
parameters.nsteps=50;
parameters.spins={'E'};
parameters.ex_prob=(operator(spin_system,{'L+'},{1})+...
                    operator(spin_system,{'L-'},{1}))/2;
parameters.ex_pump=(operator(spin_system,{'L+'},{2})+...
                    operator(spin_system,{'L-'},{2}))/2;
parameters.output='brief';
parameters.grid='rep_2ang_1600pts_sph';

% Build the time axis
time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1);

% Simulation (numerical)
deer_num=powder(spin_system,@deer_3p_hard_deer,parameters,'deer');

% Simulation (analytical)
D=xyz2dd(inter.coordinates{1},inter.coordinates{2},sys.isotopes{1},sys.isotopes{2});
D=mean(inter.zeeman.eigs{1})*mean(inter.zeeman.eigs{2})*D/(spin_system.tols.freeg^2);
deer_anl=0.35*deer_analyt(D,0,time_axis);

% Plotting (numerical)
figure(); subplot(1,2,1); 
plot(1e6*time_axis,imag(deer_num.deer_trace)); 
axis([0 0.5 -0.1 0.4]); grid ; 
xlabel('time, microseconds'); 
title('numerical result');

% Plotting (analytical)
subplot(1,2,2); plot(1e6*time_axis,deer_anl); 
axis([0 0.5 -0.1 0.4]); kgrid;
xlabel('time, microseconds'); 
title('analytical result'); 

end

