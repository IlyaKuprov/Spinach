% CLIP-HSQC spectrum of a four-spin system in a liquid crystal 
% with a user-specified order matrix.
%
% Calculation times: seconds.
%
% i.kuprov@soton.ac.uk

function rdc_fourspin()

% Magnet field
sys.magnet=5.9;

% Spin system and interactions
sys.isotopes={'1H','1H','13C','13C'};
inter.zeeman.scalar={1,10,80,50};
inter.coupling.scalar={0     2.5   250   50;
                       2.5   0     50    0;
                       250   50    0     0; 
                       50    0     0     0};
inter.coordinates={[1.0 0.0 0.0];
                   [2.0 0.0 0.0];
                   [3.0 0.0 0.0];
                   [4.0 0.0 0.0]};
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
parameters.sweep=[5000 5000];
parameters.offset=[4250 1200];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.J=140;
parameters.needs={'rdc'};
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@clip_hsqc,parameters,'nmr');

% Apodization
fid.pos=apodization(fid.pos,'sqcosbell-2d');
fid.neg=apodization(fid.neg,'sqcosbell-2d');

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 1.0 0.05 1.0],2,256,6,'negative');

end

