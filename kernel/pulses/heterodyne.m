% Filters out high frequency signals for pulse waveform as recorded by an
% oscilloscope. Supports GPU usage. Syntax:
% 
%                heterodyne(time_grid,exp_data,carrier_freq)
%
% Parameters:
%
%   time_grid      -  Time grid supplied by RLC circuit data,
%                     a column vector.
%
%   exp_data       -  Amplitude data, supplied as a column vector.
%
%   carrier_freq   -  carrier frequency, Hz
%
% Outputs:
%
%       X        - in-phase part of the rotating frame
%                  pulse waveform distorted by the RLC
%                  response, a column vector of real 
%                  numbers
%
%       Y        - out-of-phase part of the rotating
%                  frame pulse waveform distorted by
%                  the RLC response, a column vector
%                  of real numbers

function [X,Y]=heterodyne(time_grid,exp_data,carrier_freq)

% Needs a grumbler

% Carrier frequency to rad/s
omega=2*pi*carrier_freq;

% Move inputs to GPU
if canUseGPU()
    exp_data=gpuArray(exp_data); time_grid=gpuArray(time_grid);
end

% Mix with carrier frequency
X=2*exp_data.*cos(omega*time_grid);
Y=2*exp_data.*sin(omega*time_grid);
clear('exp_data','time_grid');

% Define a lowpass filter
d=designfilt('lowpassfir','SampleRate',          2.5e9, ...
                          'PassbandFrequency',   1e6,   ...
                          'StopbandFrequency',   1.5e7, ...
                          'PassbandRipple',      0.5,   ...
                          'StopbandAttenuation', 50);
B = d.Coefficients; if canUseGPU(), B=gpuArray(B); end

% Apply the lowpass filter
X=gather(fftfilt(B,X)); Y=gather(fftfilt(B,Y));

end
