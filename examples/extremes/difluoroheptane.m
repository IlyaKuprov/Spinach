% 19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by
% explicit time-domain evolution in Liouville space.
%
% WARNING: needs 32 CPU cores, 128 GB of RAM and
%          a Titan V or later.
%
% Run time on the above: minutes
%
% ilya.kuprov@weizmann.ac.il
% b.linclau@soton.ac.uk

function difluoroheptane()

% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', ...
              '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', ...
              '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'};  

% Shifts
inter.zeeman.scalar=cell(1,23);
inter.zeeman.scalar{14}=   1.0092;
inter.zeeman.scalar{15}=   1.0092;
inter.zeeman.scalar{16}=   1.0092;
inter.zeeman.scalar{21}=   1.0092;
inter.zeeman.scalar{22}=   1.0092;
inter.zeeman.scalar{23}=   1.0092;
inter.zeeman.scalar{11}=   4.6834; 
inter.zeeman.scalar{17}=   4.6834;
inter.zeeman.scalar{10}=-184.1865;
inter.zeeman.scalar{18}=-184.1865;
inter.zeeman.scalar{8}=    1.7970;
inter.zeeman.scalar{9}=    1.7970;
inter.zeeman.scalar{13}=   1.6942;
inter.zeeman.scalar{20}=   1.6942;
inter.zeeman.scalar{19}=   1.6370;
inter.zeeman.scalar{12}=   1.6370;

% J-couplings
inter.coupling.scalar=cell(23);  
inter.coupling.scalar{19, 20}=  -14.4404; 
inter.coupling.scalar{17, 8}=     2.1802;
inter.coupling.scalar{17, 9}=    10.0564;
inter.coupling.scalar{17, 19}=    4.5308;
inter.coupling.scalar{17, 20}=    7.5739;
inter.coupling.scalar{12, 13}=  -14.4404;
inter.coupling.scalar{11, 8}=    10.0564;
inter.coupling.scalar{11, 9}=     2.1802;
inter.coupling.scalar{11, 12}=    4.5308;
inter.coupling.scalar{11, 13}=    7.5739;
inter.coupling.scalar{8,  9}=   -15.1172;
inter.coupling.scalar{10, 8}=    14.0478; 
inter.coupling.scalar{10, 9}=    38.0147; 
inter.coupling.scalar{10, 11}=   49.6307;
inter.coupling.scalar{10, 13}=   18.4317; 
inter.coupling.scalar{10, 12}=   27.7945;
inter.coupling.scalar{18, 9}=    14.0478;
inter.coupling.scalar{18, 8}=    38.0147;
inter.coupling.scalar{18, 17}=   49.6307;
inter.coupling.scalar{18, 20}=   18.4317; 
inter.coupling.scalar{18, 19}=   27.7945; 
inter.coupling.scalar{10, 18}=    1.2295;
inter.coupling.scalar{12, 14}=    7.45;
inter.coupling.scalar{12, 15}=    7.45;
inter.coupling.scalar{12, 16}=    7.45;
inter.coupling.scalar{13, 14}=    7.45;
inter.coupling.scalar{13, 15}=    7.45;
inter.coupling.scalar{13, 16}=    7.45;
inter.coupling.scalar{19, 21}=    7.45;
inter.coupling.scalar{19, 22}=    7.45;
inter.coupling.scalar{19, 23}=    7.45;
inter.coupling.scalar{20, 21}=    7.45;
inter.coupling.scalar{20, 22}=    7.45;
inter.coupling.scalar{20, 23}=    7.45;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1; bas.manual=false(3,23);
bas.manual(1,[14 15 16 12 13 10 11 8 9])=1;
bas.manual(2,[12 13 10 11 8 9 17 18 19 20])=1;
bas.manual(3,[8 9 17 18 19 20 21 22 23])=1;
bas.sym_group={'S3','S3'};
bas.sym_spins={[14 15 16],[21 22 23]};
bas.zero_quantum={'1H'};
bas.projections=1;

% Greedy parallelisation
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'19F'};
parameters.rho0=state(spin_system,'L+','19F');
parameters.coil=state(spin_system,'L+','19F');
parameters.decouple={};
parameters.offset=-86700;
parameters.sweep=300;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=real(fftshift(fft(fid,parameters.zerofill)));

% Plotting
kfigure(); plot_1d(spin_system,spectrum,parameters);

end

