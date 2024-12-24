% Rapid scan ESR spectrum of a nitroxide radical.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function rapidscan_nitroxide()

% Centre field
sys.magnet=3.5;

% Spin system properties
sys.isotopes={'14N','E'};
inter.zeeman.matrix{1}=zeros(3);
inter.zeeman.matrix{2}=[2.0104    0.0000    0.0001
                        0.0000    2.0064    0.0000
                        0.0001    0.0000    2.0021];
inter.coupling.matrix{1,1}=[];
inter.coupling.matrix{2,2}=[];
inter.coupling.matrix{1,2}=[0.6178         0    0.3161
                                 0    0.5633         0
                            0.3161         0    4.1115]*1e7;
inter.coupling.matrix{2,1}=[0.6178         0    0.3161
                                 0    0.5633         0
                            0.3161         0    4.1115]*1e7;

% Simulation parameters
bas.formalism='sphten-liouv';
bas.approximation='none';
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.equilibrium='dibari';
inter.temperature=100;
inter.tau_c={2e-11};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.mw_pwr=2*pi*1e3;
parameters.sweep=[-0.011 -0.003];
parameters.nsteps=500;
parameters.timestep=1e-8;

% Run the experiment
[x_axis,spectrum]=rapidscan(spin_system,parameters);

% Plot the result
figure(); plot(x_axis,real(spectrum)); kgrid;
kxlabel('Magnetic induction, T'); xlim tight;
kylabel('Signal intensity, a.u.');

end

