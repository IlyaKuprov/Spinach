% A minimal example of H(CA)NH pulse sequence simulation.
%
% Calculation time: seconds.
%
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function hcanh_simple()

% Magnet field
sys.magnet=14.1;

% Spin system
sys.isotopes={'15N','13C','1H','13C','1H'};
sys.labels={'N','CA','H','C','HA'};

% Interactions
inter.zeeman.scalar={110 60 8 180 4};
inter.coupling.scalar=cell(5);
inter.coupling.scalar{1,3}=92;
inter.coupling.scalar{1,2}=11;
inter.coupling.scalar{1,4}=15;
inter.coupling.scalar{1,5}=1;
inter.coupling.scalar{2,4}=55;
inter.coupling.scalar{2,3}=2;
inter.coupling.scalar{2,5}=140;
inter.coupling.scalar{3,4}=4;
inter.coupling.scalar{3,5}=8;
inter.coupling.scalar{4,5}=4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H','15N','1H'};
parameters.sweep=[5000 3000 5000];
parameters.offset=[3600 -6600 3600];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hcanh,parameters,'nmr');

% Apodisation
fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.pos_neg=apodisation(spin_system,fid.pos_neg,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.neg_pos=apodisation(spin_system,fid.neg_pos,{{'sqcos'},{'sqcos'},{'sqcos'}});
fid.neg_neg=apodisation(spin_system,fid.neg_neg,{{'sqcos'},{'sqcos'},{'sqcos'}});

% F3 Fourier transform
f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3);
f3_pos_neg=fftshift(fft(fid.pos_neg,parameters.zerofill(3),3),3);
f3_neg_pos=fftshift(fft(fid.neg_pos,parameters.zerofill(3),3),3);
f3_neg_neg=fftshift(fft(fid.neg_neg,parameters.zerofill(3),3),3);

% Absorption part of F3 signal
f3_pos=f3_pos_pos+conj(f3_neg_neg);
f3_neg=f3_neg_pos+conj(f3_pos_neg);

% F2 Fourier transform
f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2);
f3f2_neg=fftshift(fft(f3_neg,parameters.zerofill(2),2),2);

% Absorption part of F2 signal
f3f2=f3f2_pos+conj(f3f2_neg);

% F1 Fourier transform
spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1);

% Plotting
kfigure(); plot_3d(spin_system,real(spectrum),parameters,...
                  10,[0.05 0.5 0.05 0.5],2,'positive');

end

