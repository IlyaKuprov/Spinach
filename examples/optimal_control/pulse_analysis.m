% An example of spectrogram analysis for a quadratic chirp 
% pulse; adapted from Matlab example set.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

function pulse_analysis()

% Quadratic chirp superposition
fs=1000; t=0:1/fs:2-1/fs;
y=chirp(t,100,1,200,'quadratic')+...
  chirp(t,200,1,100,'quadratic');

% Do the plotting
figure(); subplot(1,2,1); plot(t,y); kgrid;
kxlabel('time, seconds'); kylabel('amplitude, a.u.');
subplot(1,2,2); scale_figure([1.75 0.75]);
spectrogram(y,100,80,100,fs,'yaxis','MinThreshold',-50);
kxlabel('time, seconds'); kylabel('frequency, Hz');

end

