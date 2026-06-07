% 3D HCCH COSY experiment on a small protein fragment.
%
% Calculation time: seconds.
%
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function hcch_cosy_simple()

% Spin system
sys.isotopes={'13C','13C','13C','1H','1H'};
sys.labels={'C','CA','CB','HA','HB'};

% Interactions
sys.magnet=14.1;
inter.zeeman.scalar={180 60 40 4 2}; 
inter.coupling.scalar{2,3}=35.0;
inter.coupling.scalar{1,3}=0.3;
inter.coupling.scalar{1,2}=55.0;
inter.coupling.scalar{3,5}=130.0; 
inter.coupling.scalar{4,5}=7.0; 
inter.coupling.scalar{3,4}=6.0;
inter.coupling.scalar{2,4}=140.0;
inter.coupling.scalar{1,4}=4.0;
inter.coupling.scalar{5,5}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Sequence parameters
parameters.J_ch=140;
parameters.J_cc=35;
parameters.delta=1.1e-3;
parameters.sweep=[4000 7500 4000];
parameters.spins={'1H','13C','1H'};
parameters.offset=[2000 7000 2000];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.decouple_f3={'13C'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@hcch_cosy,parameters,'nmr');

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
kfigure(); plot_3d(spin_system,imag(spectrum),parameters,...
                   10,[0.05 0.25 0.05 0.25],2,'positive');

end

