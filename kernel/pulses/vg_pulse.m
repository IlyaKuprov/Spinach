% Veshtort-Griffin shaped pulses, generated from tables given in
%
%            http://dx.doi.org/10.1002/cphc.200400018
%
% There are good reasons to believe (see Section 2.2 of the paper) that
% these are the best possible pulses within their design specifications
% and basis sets. Syntax:
%
%             waveform=vg_pulse(pulse_name,npoints,duration)
%
% Parameters:
%
%     pulse_name - a character string, one of the following: E0A, 
%                  E0B, E100A, E100B, E200A, E200D, E200F, E300C,
%                  E300F, E400B, E300A, E500A, E500B, E500C, E600A,
%                  E600C, E600F, E800A, E800B, E1000B
%
%     npoints    - number of discrete time intervals in the pulse
%
%     duration   - duration of the pulse, seconds
%
% Outputs:
%
%     waveform   - amplitude of the pulse at each interval (there
%                  is no phase modulation), normalised to produce
%                  a 90-degree pulse, rad/s
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=vg_pulse.m>

function waveform=vg_pulse(pulse_name,npoints,duration)

% Check consistency
grumble(pulse_name,npoints,duration);

% Stockpile names
names={'E0A'   'E0B'   'E100A' 'E100B' 'E200A' 'E200D' 'E200F' 'E300C' 'E300F' 'E400B' ...
       'E300A' 'E500A' 'E500B' 'E500C' 'E600A' 'E600C' 'E600F' 'E800A' 'E800B' 'E1000B'};

% Stockpile cosine coefficients
cosines=[...
 -0.763 -0.763 -0.734  0.265  0.282  0.276  0.251  0.228  0.269  0.231  0.222  0.240  0.242  0.222 -0.739  0.239  0.249  0.238  0.235  0.254
 -0.905  0.869 -0.650  0.974 -0.314  1.213 -0.931  0.592  1.126  0.933 -1.440  2.042  0.405  0.819  0.847  1.546  1.305 -1.816 -1.192  1.239
  2.472 -0.171  2.395 -1.502 -0.534 -1.708  0.675 -1.261 -1.566 -1.183  0.164 -1.300 -0.333  0.080 -0.077 -1.609 -1.032  0.999  1.553 -1.345
 -0.811  0.039 -1.132 -0.208  0.612  0.050 -0.018  0.475 -0.114 -0.242  2.366 -1.880 -0.332 -1.473  0.531  0.363  0.388  1.235  0.122 -0.753
 -0.088 -0.005  0.058  0.498 -0.040  0.201 -0.086  0.027  0.229  0.256 -1.646  0.904 -0.257  0.248 -0.357 -0.774 -1.735  0.152 -0.716  0.797
  0.001  0.003  0.036  0.158  0.084  0.109  0.070 -0.069  0.116  0.041  0.246  0.009  0.260  0.128 -0.140  0.078  1.087 -0.710 -0.336  0.610
  0.107  0.001  0.053 -0.165 -0.157 -0.082 -0.026  0.008 -0.015  0.025  0.192  0.072  0.064  0.018 -0.183  0.208 -0.273 -0.002  0.525 -1.169
 -0.054  0.001 -0.045 -0.079  0.002 -0.029  0.012  0.011 -0.058 -0.052 -0.211  0.032 -0.014  0.044  0.087  0.071  0.115 -0.120 -0.185  0.377
  0.006  0.001  0.011  0.051  0.072  0.012  0.000 -0.005  0.006  0.002  0.118 -0.110 -0.014 -0.029  0.090 -0.057 -0.003  0.011  0.022  0.103
  0.002  0.001  0.000  0.035 -0.021  0.009  0.001 -0.002  0.010  0.007  0.003  0.033 -0.034 -0.024  0.012 -0.009 -0.036 -0.032  0.008 -0.055
  0.004  0.000  0.000 -0.015 -0.003 -0.005  0.002  0.001  0.003  0.001 -0.058  0.006  0.021  0.007 -0.029 -0.002 -0.017  0.083 -0.041  0.009
 -0.002  0.000  0.000 -0.015 -0.005 -0.004  0.000 -0.001 -0.004 -0.001  0.043 -0.009 -0.004  0.001 -0.021 -0.021  0.026  0.064  0.045 -0.043
  0.001  0.000 -0.001  0.003  0.006  0.000  0.001  0.000 -0.002 -0.002 -0.008  0.008  0.005  0.000 -0.003  0.015 -0.026 -0.051 -0.007  0.020
  0.001  0.000  0.001  0.006  0.002  0.000  0.001  0.000  0.000  0.000 -0.008 -0.007 -0.003  0.000  0.008 -0.001  0.020 -0.038 -0.021  0.013
  0.001  0.000  0.000 -0.001 -0.002 -0.001  0.000  0.000  0.000  0.000  0.010 -0.001 -0.002 -0.002  0.005 -0.001 -0.012  0.011  0.013 -0.008
  0.000  0.000  0.000 -0.002  0.000 -0.001  0.000  0.000  0.000  0.000 -0.004  0.001  0.001 -0.001 -0.002 -0.002  0.001  0.008  0.001 -0.020
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001 -0.002 -0.001  0.000 -0.003  0.000  0.000  0.000 -0.005  0.020
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.003  0.000  0.001 -0.001 -0.001 -0.002 -0.003 -0.001  0.003 -0.007
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001 -0.001 -0.001  0.000  0.000  0.000  0.001  0.001 -0.002 -0.001
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001  0.000  0.000  0.000  0.000 -0.001 -0.001  0.000  0.000
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.001  0.000  0.000  0.000  0.000  0.000 -0.001 -0.002  0.001  0.000];

% Stockpile sine coefficients
sines=[...
  0.000  0.000  0.092 -2.104  1.459 -1.270 -1.273 -1.279 -0.788 -0.533 -0.772 -0.681 -3.105 -0.215 -0.360  0.249 -2.489 -0.238 -0.773  0.242
  0.000  0.000 -0.729 -0.025  1.934 -0.244  0.913  0.828 -1.010 -1.288 -0.064  0.389  0.674 -1.555 -3.153 -1.861  0.874 -0.228  1.045 -1.607
  0.000  0.000  0.564  0.795 -2.256  0.615 -0.312  0.029  0.875  1.006  0.427 -0.552 -1.212  0.357  0.210 -0.773 -0.847 -2.724 -0.918 -2.102
  0.000  0.000 -0.057  0.222  0.189  0.106  0.106 -0.082  0.078  0.022 -0.269  0.177  1.609  0.525  1.867  1.583  0.610  0.724 -1.052  1.532
  0.000  0.000 -0.032 -0.292  0.234 -0.125 -0.023 -0.038 -0.033 -0.010  0.284  0.313 -0.346 -0.034 -0.113 -0.285  0.255  1.886  1.391 -0.409
  0.000  0.000 -0.063 -0.113  0.037 -0.056 -0.003  0.047 -0.102 -0.061 -0.144 -0.193  0.025  0.013 -0.168  0.121 -0.123 -0.316 -0.271  1.110
  0.000  0.000  0.079  0.093 -0.025  0.036  0.009 -0.012  0.015 -0.007 -0.064  0.043 -0.017 -0.029 -0.081 -0.039  0.010 -0.305 -0.101 -0.359
  0.000  0.000 -0.033  0.053 -0.063  0.017 -0.005 -0.002  0.024  0.026  0.110 -0.020 -0.085 -0.039 -0.091 -0.017 -0.003 -0.040  0.076 -0.091
  0.000  0.000  0.004 -0.028  0.036 -0.006  0.003  0.003  0.001 -0.001 -0.047 -0.017  0.085  0.018  0.050 -0.068 -0.077 -0.020 -0.035  0.071
  0.000  0.000 -0.001 -0.024  0.007 -0.006 -0.001  0.001 -0.007 -0.001 -0.008  0.035 -0.034  0.009  0.066  0.054  0.090  0.026  0.018  0.056
  0.000  0.000  0.004  0.007 -0.010  0.002  0.000  0.000 -0.003 -0.002  0.028 -0.013  0.022 -0.002 -0.015 -0.005 -0.053  0.008  0.027 -0.102
  0.000  0.000 -0.003  0.009  0.000  0.002  0.000  0.001  0.001  0.001 -0.024 -0.002 -0.012  0.000 -0.018  0.001  0.033  0.022 -0.042  0.025
  0.000  0.000  0.001 -0.002  0.001  0.000  0.000  0.001  0.000  0.001  0.008  0.003 -0.001 -0.002 -0.003  0.001 -0.014  0.013  0.016  0.027
  0.000  0.000  0.000 -0.004  0.001  0.000  0.000  0.000 -0.001  0.000  0.004 -0.003  0.002 -0.001  0.001 -0.002  0.002 -0.017  0.008 -0.010
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001  0.000 -0.007  0.002 -0.002  0.000  0.004 -0.002  0.000 -0.008 -0.010  0.003
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.004  0.000  0.001  0.000  0.001  0.001  0.001  0.002  0.003 -0.005
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001 -0.002 -0.001 -0.002  0.000 -0.002  0.002  0.001  0.000
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001  0.000  0.000 -0.001 -0.001 -0.001  0.002  0.003 -0.002  0.004
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.001  0.000 -0.001 -0.001  0.000  0.000 -0.002  0.000  0.001 -0.002
  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 -0.001  0.000 -0.001 -0.001  0.000 -0.001  0.001 -0.001 -0.001 -0.002];

% Parse the pulse name
pulse_number=find(strcmp(names,pulse_name));
if isempty(pulse_number), error('incorrect pulse name.'); end

% Preallocate the pulse
waveform=zeros(1,npoints);

% Build the time grid
time_grid=linspace(0,2*pi,npoints);

% Compute cosine terms
for k=0:20
    waveform=waveform+cosines(k+1,pulse_number)*cos(k*time_grid);
end

% Compute sine terms
for k=1:20
    waveform=waveform+sines(k,pulse_number)*sin(k*time_grid);
end

% Normalise the waveform
waveform=2*pi*waveform/duration;

end

% Consistency enforcement
function grumble(pulse_name,npoints,duration)
if ~ischar(pulse_name)
    error('pulse_name must be a character string.');
end
if (numel(npoints)~=1)||(~isnumeric(npoints))||(~isreal(npoints))||...
   (npoints<1)||(mod(npoints,1)~=0)
    error('npoints must be a positive real integer greater than 1.');
end
if (numel(duration)~=1)||(~isnumeric(duration))||...
   (~isreal(duration))||(duration<0)
    error('duration must be a positive real number.');
end
end

% It is not the works, but the belief which is here decisive and deter-
% mines the order of rank - to employ once more an old religious formu-
% la with a new and deeper meaning - it is some fundamental certainty 
% which a noble soul has about itself, something which is not to be so-
% ught, is not to be found, and perhaps, also, is not to be lost. The 
% noble soul has reverence for itself.
%
% Friedrich Nietzsche, in "Beyond Good and Evil"

