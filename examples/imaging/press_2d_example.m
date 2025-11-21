% 2D PRESS example. Three independent spin systems are localised 
% in three spots of a 2D sample. The spots are slectively excited
% and their NMR spectra recorded. The followng are the frequencies
% to excite the three substances:
%
%    parameters.rf_frq_list={-120e3 -100e3}  - substance A
%
%    parameters.rf_frq_list={-80e3 -10e3}    - substance B
%
%    parameters.rf_frq_list={+30e3 +100e3}   - substance C
%
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function press_2d_example()

% Magnetic induction
sys.magnet=3.0;

% Spin systems
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
parameters.sweep=1000;
parameters.axis_units='ppm';
parameters.npoints=128;
parameters.invert_axis=1;
parameters.ss_grad_amp=[25e-3 25e-3];
parameters.image_size=[101 105];    
parameters.sp_grad_amp=5e-3; % T/m
parameters.sp_grad_dur=5e-4;

% Pulse parameters
parameters.rf_phi={pi/2 pi/2};
parameters.rf_frq_list={-120e3 -100e3};       % change this to scan through the sample
parameters.rf_amp_list={2*pi*5000 2*pi*5000};
parameters.rf_dur_list={0.5e-4 1.0e-4};
parameters.max_rank={3 3};

% Sample geometry
parameters.dims=[0.30 0.25];
parameters.npts=[108 90];
parameters.deriv={'period',3};

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
[X,Y]=ndgrid(1:parameters.npts(1),1:parameters.npts(2));
parameters.rho0_ph={(X-15).^2+(Y-15).^2<10^2;
                    (X-30).^2+(Y-40).^2<10^2;
                    (X-70).^2+(Y-80).^2<10^2};   
parameters.rho0_st={state(spin_system,'Lz',[1 2]);                             
                    state(spin_system,'Lz',[3 4]);
                    state(spin_system,'Lz',[5 6])};
parameters.coil_ph={ones(parameters.npts)};
parameters.coil_st={state(spin_system,'L+','all')};

% Show the phantom
kfigure(); scale_figure([2.0 1.0]); 
subplot(1,3,1); 
mri_2d_plot(parameters.rho0_ph{1}+...
            parameters.rho0_ph{2}+...
            parameters.rho0_ph{3},parameters,'phantom');
ktitle('sample phantom'); drawnow();

% Run active volume diagnostics
phan=imaging(spin_system,@press_voxel_2d,parameters);
subplot(1,3,2); mri_2d_plot(phan,parameters,'phantom');
ktitle('active volume'); drawnow();

% Run PRESS 2D
fid=imaging(spin_system,@press_2d,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'}});

% Magnitude Fourier transform
spec=abs(fftshift(fft(ifftshift(fid))));

% Plotting
subplot(1,3,3); plot_1d(spin_system,spec,parameters);
ktitle('active volume spectrum');

end

