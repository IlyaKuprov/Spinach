% Field swept EPR spectrum of nitroxide, computed by finding
% resonance fields and transition moments.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

function fieldsweep_nitroxide()

% Isotopes                          
sys.isotopes={'E','14N'};
                          
% Magnet field (must be 1)
sys.magnet=1;

% Interactions
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=[2.01045  0.00000  0.00000
                        0.00000  2.00641  0.00000
                        0.00000  0.00000  2.00211];
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=[1.2356  0.0000  0.6322
                            0.0000  1.1266  0.0000
                            0.6322  0.0000  8.2230]*1e7;

% Temperature
inter.temperature=298;
                        
% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.grid=6;
parameters.mw_freq=9e9;
parameters.fwhm=1e-5;
parameters.int_tol=1e-5;
parameters.tm_tol=0.1;
parameters.window=[0.316 0.326];
parameters.npoints=1024;
parameters.rspt_order=Inf;

% Run the simulation
[b_axis,spec]=fieldsweep(spin_system,parameters);

% Plotting
figure(); plot(b_axis',spec');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
axis tight; kgrid;

end

