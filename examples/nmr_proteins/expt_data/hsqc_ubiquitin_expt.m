% Experimental HSQC spectrum of human ubiquitin.
%
% Donghan Lee (Max Planck Institute)
% Ilya Kuprov (University of Southampton)

function hsqc_ubiquitin_expt()

% Load the file
load('hsqc_ubiquitin_expt.mat','fid');
spin_system.inter.magnet=11.7395;
parameters.zerofill=[1024 1024];
parameters.sweep=[2000 4000];
parameters.offset=[-5870 3753];
parameters.spins={'15N','1H'};
parameters.axis_units='ppm';

% Phasing
phase_f1=exp(-1i*0.7);

% Apodisation
fid.pos=apodisation(spin_system,fid.pos,{{'cos'},{'cos'}})*phase_f1;
fid.neg=apodisation(spin_system,fid.neg,{{'cos'},{'cos'}})*phase_f1;

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form States signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
spectrum=flip(flip(spectrum,1),2);
spin_system.sys.disable={'colorbar'}; 
spin_system.sys.output=1; 
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-imag(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'positive');

end

