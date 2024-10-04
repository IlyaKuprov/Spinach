% Phosphorus system simulation for Gerhard Hagele. Done by brute
% force Liouville space time propagation.
%
% WARNING: needs 32 CPU cores, 128 GB of RAM and
%          a Titan V or later.
%
% Run time on the above: hours
%
% i.kuprov@soton.ac.uk

function phosphorus_cluster()

% Magnet induction 
sys.magnet=7.04;

% Isotopes
sys.isotopes={'31P','31P','31P','31P','31P','31P','31P',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar=cell(1,34);
inter.zeeman.scalar{1}=-99.6;
inter.zeeman.scalar{2}=-0.33;
inter.zeeman.scalar{3}=-0.33;
inter.zeeman.scalar{4}=-0.33;
inter.zeeman.scalar{5}=-156.7;
inter.zeeman.scalar{6}=-156.7;
inter.zeeman.scalar{7}=-156.7;
for n=8:34
    inter.zeeman.scalar{n}=0.21;
end

% J-couplings
inter.coupling.scalar=cell(34,34);  
inter.coupling.scalar{1,2}=-323.22;
inter.coupling.scalar{1,3}=-323.22;
inter.coupling.scalar{1,4}=-323.22;
inter.coupling.scalar{1,5}=  46.18;
inter.coupling.scalar{1,6}=  46.18;
inter.coupling.scalar{1,7}=  46.18;
inter.coupling.scalar{2,5}=-354.82;
inter.coupling.scalar{3,6}=-354.82;
inter.coupling.scalar{4,7}=-354.82;
inter.coupling.scalar{2,6}=  25.79;
inter.coupling.scalar{3,7}=  25.79;
inter.coupling.scalar{4,5}=  25.79;
inter.coupling.scalar{2,7}=  -9.04;
inter.coupling.scalar{3,5}=  -9.04;
inter.coupling.scalar{4,6}=  -9.04;
inter.coupling.scalar{2,3}= -16.63;
inter.coupling.scalar{3,4}= -16.63;
inter.coupling.scalar{4,2}= -16.63;
inter.coupling.scalar{5,6}=-214.10;
inter.coupling.scalar{6,7}=-214.10;
inter.coupling.scalar{7,5}=-214.10;
for n=8:16
    inter.coupling.scalar{2,n}=4.0;
end
for n=17:25
    inter.coupling.scalar{3,n}=4.0;
end
for n=26:34
    inter.coupling.scalar{4,n}=4.0;
end

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;
bas.longitudinals={'1H'};
bas.projections=1;

% Symmetry
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[8 9 10],[17 18 19],[26 27 28]};

% Use greedy parallelisation
sys.enable={'greedy','gpu'};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'31P'};
parameters.rho0=state(spin_system,'L+','31P');
parameters.coil=state(spin_system,'L+','31P');
parameters.decouple={};
parameters.offset=-10000;
parameters.sweep=30000;
parameters.npoints=16536;
parameters.zerofill=65536;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',5}});

% Fourier transform
spectrum=real(fftshift(fft(fid,parameters.zerofill)));

% Plotting
figure(); plot_1d(spin_system,spectrum,parameters);

end

