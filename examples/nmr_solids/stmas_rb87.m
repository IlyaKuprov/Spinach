% Shifted-echo 87Rb satellite-transition MAS with finite excitation,
% conversion, and central-transition inversion pulses.
% Pulse schedule and sample calibration: DOI 10.1016/j.jmro.2026.100218.
%
% talos@spindynamics.org

function stmas_rb87()

% System and quadrupolar coupling
sys.magnet=2*pi*304.47e6/spin('1H');
sys.isotopes={'87Rb'};
sys.parallel={'processes',64};
inter.zeeman.scalar={4700e6*2*pi/(sys.magnet*spin('87Rb'))};
inter.zeeman.matrix={diag([-250 -250 500])};
inter.coupling.matrix={eeqq2nqi(4e6,0,3/2,[0 0 0])};

% Basis set
bas.formalism='zeeman-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse sequence and rotor-synchronous sampling
parameters.rate=50e3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=15;
parameters.grid='rep_2ang_400pts_sph';
parameters.npoints=[256 512];
parameters.tau=parameters.npoints(2)/(2*parameters.rate);
parameters.spins={'87Rb'};
parameters.offset=0;
parameters.rframes={{'87Rb',2}};
parameters.rho0=state(spin_system,'Lz',1);
parameters.coil=coil_state(spin_system,'CT-',1);
parameters.pulse_amp=2*pi*[100e3 100e3 10e3];
parameters.pulse_dur=[1.6 1.6 25]*1e-6;

% Simulate the complex echo in direct-by-indirect sample order
fid=singlerot(spin_system,@stmas,parameters,'labframe');

% Process both time dimensions with their distinct spectral widths
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}});
spectrum=fftshift(fft2(fid));
parameters.sweep=parameters.rate*[17/9 1];
parameters.offset=[0 0];
parameters.axis_units='kHz';

% Plot the phase-insensitive shifted-echo spectrum
kfigure();
plot_2d(spin_system,abs(spectrum),parameters,20,[0.05 0.5 0.05 0.5],2,256,6,'positive');

end


