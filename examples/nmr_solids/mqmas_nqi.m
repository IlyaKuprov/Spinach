% Rotor-synchronous MQMAS spectrum of a 87Rb compound,
% transmitter set to the isotropic chemical shift.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function mqmas_nqi()

% System specification: just the NQI
sys.magnet=9.4; sys.isotopes={'87Rb'};
inter.coupling.matrix{1,1}=eeqq2nqi(5e6,0.50,3/2,[0 0 0]);

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=62.5e3;
parameters.axis=[1 1 1];
parameters.max_rank=7;
parameters.grid='rep_2ang_1600pts_sph';
parameters.npoints=[128 128];
parameters.zerofill=[256 256];
parameters.sweep=parameters.rate;
parameters.offset=0;
parameters.spins={'87Rb'};
parameters.rframes={{'87Rb',2}};
parameters.mq_order=3;
parameters.axis_units='ppm';
parameters.rho0=state(spin_system,'Lz','87Rb');
parameters.coil=state(spin_system,'L+','87Rb');
parameters.pulse_amp=2*pi*[250e3 250e3];
parameters.pulse_dur=[2e-6 1e-6];

% Simulation
fid=singlerot(spin_system,@mqmas,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}});
    
% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.05 0.5 0.05 0.5],2,256,6,'positive');
                                     
end

