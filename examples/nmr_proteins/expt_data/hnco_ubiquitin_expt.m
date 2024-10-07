% Experimental HNCO spectrum of human ubiquitin.
%
% Donghan Lee (Max Planck Institute)
% Ilya Kuprov (University of Southampton)

function hnco_ubiquitin_expt()

% Data loading and truncation
load('hnco_ubiquitin_expt.mat','fid'); 
fid=fid(1:64,:,:); 

% Apodisation
fid=apodisation([],fid,{{'cos'},{'cos'},{'cos'}});

% F3 processing
fid=fft(fid,256,1); fid=fid*exp(-1i*0.75);

% F2 processing
fid=real(fid(:,1:2:end,:))+1i*real(fid(:,2:2:end,:));
fid=fft(fid,256,2); 

% F1 processing
fid=real(fid(:,:,1:2:end))-1i*real(fid(:,:,2:2:end));
spectrum=fft(fid,256,3);

% Window shifting
spectrum=permute(spectrum,[3 2 1]);
spectrum=fftshift(spectrum,1);
spectrum=fftshift(spectrum,2);
spectrum=fftshift(spectrum,3);

% Baseline correction
spectrum=spectrum-spectrum(end,end,end);

% Water signal elimination
spectrum(:,:,1:40)=0;

% Parameters
spin_system.inter.magnet=11.7395;
spin_system.sys.disable={'colorbar'};
parameters.spins={'15N','13C','1H'};
parameters.offset=[-5870 22164 3653];
parameters.sweep=[2000 1500 3000];
parameters.zerofill=[256 256 256];
parameters.npoints=[64 64 64];
parameters.axis_units='ppm';

% Plotting
spin_system.sys.disable={'colorbar'}; spin_system.sys.output=1; figure();
plot_3d(spin_system,real(spectrum),parameters,10,[0.2 0.9 0.2 0.9],2,'positive');

end

