% A minimal HNCOCA pulse sequence simulation.
%
% Calculation time: seconds.
%
% ledwards@cbs.mpg.de
% ilya.kuprov@weizmann.ac.il

function hncoca_simple()

% Magnet field
sys.magnet=14.1;

% Spin system
sys.isotopes={'15N','13C','1H','13C'};
sys.labels={'N','CA','H','C'};

% Interactions
inter.zeeman.scalar={110 55 8 180};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,3}=92;
inter.coupling.scalar{1,2}=11;
inter.coupling.scalar{1,4}=15;
inter.coupling.scalar{2,4}=55;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tau=[2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3];
parameters.spins={'15N','13C','1H'};
parameters.offset=[-7200 5600 4800];
parameters.sweep=[5000 8000 5000];
parameters.npoints=[63 64 65];
parameters.zerofill=[255 256 257];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hncoca,parameters,'nmr');

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
                  10,[0.2 0.9 0.2 0.9],2,'positive');

end

