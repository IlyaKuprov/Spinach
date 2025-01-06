% Chirped inversion pulse.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function shaped_pulse_chirp_xy()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Zeeman interactions
inter.zeeman.scalar=num2cell(linspace(-4,4,31));

% Couplings
inter.coupling.scalar=cell(31);
for n=1:30
    inter.coupling.scalar{n,n+1}=20;
end

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'Lz','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=5100;
parameters.npoints=2048;
parameters.zerofill=8192;
parameters.axis_units='Hz';

% Pulse infrastructure
H=hamiltonian(assume(spin_system,'nmr'));
R=relaxation(spin_system);
K=kinetics(spin_system);
Lp=operator(spin_system,'L+','1H');
Lm=operator(spin_system,'L-','1H');
Lx=(Lp+Lm)/2; Ly=(Lp-Lm)/2i;

% Chirp waveform in amplitude-frequency coordinates
[Cx,Cy,durs]=chirp_pulse(500,0.1,2000,20,'wurst-adaptive');

% Soft pulse
parameters.rho0=shaped_pulse_xy(spin_system,H,{Lx,Ly},{Cx,Cy},...
                                durs,parameters.rho0,'expv-pwc');
                            
% Homospoil
parameters.rho0=homospoil(spin_system,parameters.rho0,'destroy');

% Global hard pulse
parameters.rho0=step(spin_system,Ly,parameters.rho0,pi/2);

% Acquisition
fid=acquire(spin_system,parameters,H,R,K);

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));
spectrum=real(spectrum);

% Plotting
figure(); scale_figure([2.0 1.0]);
subplot(1,2,1); plot(cumsum(durs),[Cx; Cy]); 
kxlabel('time, seconds'); ktitle('chirp pulse'); 
kylabel('amplitude, rad/s'); kgrid; xlim tight;
subplot(1,2,2); plot_1d(spin_system,spectrum,parameters);
ktitle('band-selective inversion');

end

