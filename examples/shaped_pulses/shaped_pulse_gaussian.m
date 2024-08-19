% Gaussian 90-degree pulse on a chain of 31 strongly coupled protons.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function shaped_pulse_gaussian()

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

% Assumptions
spin_system=assume(spin_system,'nmr');

% Hamiltonian superoperator
H=hamiltonian(spin_system);

% Control and offset operators
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');

% Initial state
rho=state(spin_system,'Lz','1H');

% Gaussian pulse setup
npoints=80; duration=0.015; 
time_grid=duration*ones(1,npoints)/npoints;
[amplitudes,phases]=read_wave('gaussian_1000.pk',npoints);

% Pulse calibration
integral=sum(amplitudes);
calibration_factor=(pi/2)*npoints/(integral*duration);
amplitudes=calibration_factor*amplitudes;

% Pulse transformation from (A,phi) to (X,Y)
[Cx,Cy]=polar2cartesian(amplitudes,phases);

% Pulse execution with a 480 Hz frequency offset
rho=shaped_pulse_xy(spin_system,H+2*pi*480*Lz,...
                    {Lx,Ly},{Cx,Cy},time_grid,rho,'expv-pwc');

% Set up acquisition
parameters.spins={'1H'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=5000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='Hz';

% Run acquisition
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,imag(spectrum),parameters);

end

