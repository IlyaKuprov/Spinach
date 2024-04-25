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
%   carrier_freq   -  carrier frequency of the instrument
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

% Get carrier frequency in rad/s
omega=2*pi*carrier_freq;

% Mix with carrier frequency on GPU if possible 
if canUseGPU()
    X=gather(2*gpuArray(exp_data).*cos(omega*gpuArray(time_grid)));
    Y=gather(2*gpuArray(exp_data).*sin(omega*gpuArray(time_grid)));
else
    X=2*exp_data.*cos(omega*time_grid);
    Y=2*exp_data.*sin(omega*time_grid);
end

% Define filter
d=designfilt('lowpassfir','SampleRate',2.5e9, ...
             'PassbandFrequency',1e6,'StopbandFrequency',1.5e7, ...
             'PassbandRipple',0.5,'StopbandAttenuation',50);
B = d.Coefficients;

% Apply filter on GPU if available
if canUseGPU()
    X=gather(fftfilt(gpuArray(B),gpuArray(X)));
    Y=gather(fftfilt(gpuArray(B),gpuArray(Y)));
else 

% Apply filter on CPU
X=fftfilt(B,X);
Y=fftfilt(B,Y);

end 

end
