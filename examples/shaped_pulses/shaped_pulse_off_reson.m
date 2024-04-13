% An off-resonance rectangular soft pulse simulated using
% the Fokker-Planck formalism.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function shaped_pulse_off_reson()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Zeeman interactions
inter.zeeman.scalar=num2cell(linspace(-4,4,31));

% Couplings
inter.coupling.scalar=cell(31);
for n=1:30
    inter.coupling.scalar{n,n+1}=10;
end

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse infrastructure
H=hamiltonian(assume(spin_system,'nmr'));
Lp=operator(spin_system,'L+','1H');
Lm=operator(spin_system,'L-','1H');
Lx=(Lp+Lm)/2; Ly=(Lp-Lm)/2i;
rho=state(spin_system,'Lz','1H');

% Soft pulse
rho=shaped_pulse_af(spin_system,H,Lx,Ly,rho,1922.4,50.0,5e-3,-pi/2,2);

% Set up acquisition
parameters.spins={'1H'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=7000;
parameters.npoints=2048;
parameters.zerofill=8192;
parameters.axis_units='Hz';

% Run acquisition
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Phasing
fid=fid*exp(-1i*0.67);

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

