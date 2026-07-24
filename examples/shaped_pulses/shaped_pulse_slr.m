% Shinnar-Le Roux band-selective 90-degree pulse on a chain of
% 31 strongly coupled protons.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function shaped_pulse_slr()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes=repmat({'1H'},1,31);

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
spin_system=assume(spin_system,'nmr');

% Hamiltonian and control operators
H=hamiltonian(spin_system);
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');

% Initial state
rho=state(spin_system,'Lz','1H');

% SLR pulse waveform
[Cx,Cy,durs]=slr_pulse(128,15e-3,4,pi/2,0.01,0.01);

% Pulse execution
rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{Cx,Cy},...
                    durs,rho,'expv-pwc');

% Set up acquisition
parameters.spins={'1H'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=5000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='Hz';

% Run acquisition
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); scale_figure([2.0 1.0]);
subplot(1,2,1); plot(cumsum(durs),[Cx;Cy]);
kxlabel('time, seconds'); ktitle('SLR pulse');
kylabel('amplitude, rad/s'); kgrid; xlim tight;
subplot(1,2,2); plot_1d(spin_system,imag(spectrum),parameters);
ktitle('band-selective excitation');

end


