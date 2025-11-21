% 1D PRESS example. Three independent spin systems are localised 
% in three areas of a 1D sample. The areas are slectively excited
% and their NMR spectra recorded. The followng are the frequencies
% to excite the three substances:
% 
%  pulse_frq=+100e3  - substances B and C
%  
%  pulse_frq=0;      - substance C
%
%  pulse_frq=-100e3  - substances A and C
%
% Simulation time: seconds, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function press_1d_example()

% Magnetic induction
sys.magnet=3.0;
sys.isotopes={'1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={-3,3,-2,2,-1,1};
inter.coupling.scalar=cell(6);
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{3,4}=20;
inter.coupling.scalar{5,6}=30;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Disable path tracing
sys.disable={'pt'};

% This needs a GPU
% sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0;
parameters.sweep=500000;
parameters.npoints=128;
parameters.axis_units='kHz';
parameters.invert_axis=1;
parameters.ss_grad_amp=30e-3;
parameters.ro_grad_amp=30e-3;

% Pulse parameters
parameters.rf_phi=pi/2;
parameters.rf_frq_list=-100e3;       % change this to scan through the sample
parameters.rf_amp_list=2*pi*5000;
parameters.rf_dur_list=0.5e-4;
parameters.max_rank=3;

% Sample geometry
parameters.dims=0.30;
parameters.npts=100;
parameters.deriv={'period',3};

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts,1)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={[zeros(10,1); ones(30,1); zeros(60,1)];
                    [zeros(60,1); ones(30,1); zeros(10,1)];
                    [ones(60,1);  ones(30,1); ones(10,1)]};
parameters.rho0_st={state(spin_system,{'Lz'},{1})+state(spin_system,{'Lz'},{2});
                    state(spin_system,{'Lz'},{3})+state(spin_system,{'Lz'},{4});
                    state(spin_system,{'Lz'},{5})+state(spin_system,{'Lz'},{6})};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','all')};

% Show the phantom
kfigure(); scale_figure([1.50 0.75]); 
subplot(1,3,1); 
plot(parameters.rho0_ph{1}+...
     parameters.rho0_ph{2}+...
     parameters.rho0_ph{3}); ylim([0 3]);
ktitle('sample phantom'); kgrid;
kxlabel('pixels'); drawnow();

% Run voxel selection diagnostics
phan=imaging(spin_system,@press_voxel_1d,parameters);

% Plotting
subplot(1,3,2); 
plot(real(phan)); kgrid;
ktitle('excitation profile'); 
kxlabel('pixels'); drawnow();

% Run PRESS 1D
parameters.sweep=1000;
parameters.axis_units='ppm';
fid=imaging(spin_system,@press_1d,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'}});

% Magnitude Fourier transform
spec=abs(fftshift(fft(ifftshift(fid))));

% Plotting
subplot(1,3,3); plot_1d(spin_system,spec,parameters);
ktitle('voxel spectrum');

end

