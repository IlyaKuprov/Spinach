% CLIP-HSQC spectrum of a C-H system in a liquid crystal 
% with a user-specified order matrix.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function rdc_twospin()

% Spin system parameters
sys.magnet=5.9;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={5.0 65.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=140;
inter.coordinates={[0.0 0.0 0.0];
                   [0.6 0.7 0.8]};
inter.order_matrix={diag([1e-3 2e-3 -3e-3])};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'colorbar'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.sweep=[3000 1000];
parameters.offset=[4250 1200];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.axis_units='ppm';
parameters.decouple_f1={};
parameters.decouple_f2={'13C'};
parameters.needs={'rdc'};
parameters.J=140;

% Simulation
fid=liquid(spin_system,@clip_hsqc,parameters,'nmr');

% Apodisation
fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}});
fid.neg=apodisation(spin_system,fid.neg,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 1.0 0.05 1.0],2,256,6,'negative');

end

