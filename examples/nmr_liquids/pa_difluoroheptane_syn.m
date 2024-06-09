% Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a
% manual basis set specification as a merger of Lie algebras of 
% the user-specified structral fragments followed by symmetry fac-
% torisation and conservation law screening. See our paper:
%
%           https://doi.org/doi/10.1021/acs.joc.4c00670
%
% for further information.
%
% Calculation time: minutes, faster with a GPU.
%
% i.kuprov@soton.ac.uk

function pa_difluoroheptane_syn()

% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', ...
              '1H', '1H',  '19F', '1H', '1H', '1H', '1H', '1H', ...
              '1H', '19F', '1H',  '1H', '1H', '1H', '1H', '1H'};

% Chemical shifts
inter.zeeman.scalar=cell(1,23);
inter.zeeman.scalar{14}= 1.0189;
inter.zeeman.scalar{15}= 1.0189;
inter.zeeman.scalar{16}= 1.0189;
inter.zeeman.scalar{21}= 1.0189;
inter.zeeman.scalar{22}= 1.0189;
inter.zeeman.scalar{23}= 1.0189;
inter.zeeman.scalar{11}= 4.6138;
inter.zeeman.scalar{18}= 4.6138;
inter.zeeman.scalar{17}= 0.0000; % actually -181.9705, but does not 
inter.zeeman.scalar{10}= 0.0000; % matter here, and faster when zero
inter.zeeman.scalar{8}=  2.0878;
inter.zeeman.scalar{9}=  1.8445;
inter.zeeman.scalar{13}= 1.7132;
inter.zeeman.scalar{19}= 1.7132;
inter.zeeman.scalar{20}= 1.7018;
inter.zeeman.scalar{12}= 1.7018;

% J-couplings
inter.coupling.scalar=cell(23);  
inter.coupling.scalar{12,13}=-13.8188; 
inter.coupling.scalar{19,20}=-13.8188;
inter.coupling.scalar{9,11}=   4.8881;
inter.coupling.scalar{9,18}=   4.8881;
inter.coupling.scalar{8,11}=   7.0776;
inter.coupling.scalar{8,18}=   7.0776;
inter.coupling.scalar{11,12}=  4.4430;
inter.coupling.scalar{18,20}=  4.4430;
inter.coupling.scalar{11,13}=  7.7255;
inter.coupling.scalar{18,19}=  7.7255;
inter.coupling.scalar{8,9}=  -14.7896;
inter.coupling.scalar{8,10}=  18.3289; 
inter.coupling.scalar{8,17}=  18.3289; 
inter.coupling.scalar{9,10}=  25.8358;
inter.coupling.scalar{9,17}=  25.8358; 
inter.coupling.scalar{10,11}= 48.5693;
inter.coupling.scalar{17,18}= 48.5693;
inter.coupling.scalar{10,13}= 17.2030;
inter.coupling.scalar{17,19}= 17.2030;
inter.coupling.scalar{10,12}= 30.4926; 
inter.coupling.scalar{17,20}= 30.4926; 
inter.coupling.scalar{10,17}=  1.8459;
inter.coupling.scalar{12,14}=  7.45;
inter.coupling.scalar{12,15}=  7.45;
inter.coupling.scalar{12,16}=  7.45;
inter.coupling.scalar{13,14}=  7.45;
inter.coupling.scalar{13,15}=  7.45;
inter.coupling.scalar{13,16}=  7.45;
inter.coupling.scalar{19,21}=  7.45;
inter.coupling.scalar{19,22}=  7.45;
inter.coupling.scalar{19,23}=  7.45;
inter.coupling.scalar{20,21}=  7.45;
inter.coupling.scalar{20,22}=  7.45;
inter.coupling.scalar{20,23}=  7.45;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1; bas.manual=false(3,23);
bas.manual(1,[14 15 16 12 13 10 11 8 9])=1;
bas.manual(2,[12 13 10 11 8 9 17 18 19 20])=1;
bas.manual(3,[8 9 17 18 19 20 21 22 23])=1;
bas.sym_group={'S3','S3'};
bas.sym_spins={[14 15 16],[21 22 23]};
bas.longitudinals={'19F'};
bas.projections=1;

% GPU is useful here
sys.enable={'gpu'};

% ZTE not useful here
sys.disable={'zte'};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=1400;
parameters.sweep=2500;
parameters.npoints=4096;
parameters.zerofill=16536;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',5.0);

% Fourier transform
spectrum=real(fftshift(fft(fid,parameters.zerofill)));

% Plotting
plot_1d(spin_system,spectrum,parameters);

end

