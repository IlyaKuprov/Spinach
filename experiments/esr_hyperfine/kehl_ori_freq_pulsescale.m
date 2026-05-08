%KEHL_ORI_FREQ_PULSESCALE Microwave pulse excitation factor for frequency-domain selection.
%
%   Spinach architecture migration May 2026 Talos

function scalefactor=kehl_ori_freq_pulsescale(parameters,DeltaOm)

% Check consistency
grumble(parameters,DeltaOm);

    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % parameters: Kehl ENDOR context parameters
    % DeltaOm: offset of the actual freq from the resonance freq
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %
    data=parameters.pulse_file;
    %%
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y(1:70),(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isfield(parameters,'multipulses')&&parameters.multipulses==true
       Y=y;
    end
%%

    %Delta omega in MHz
    DeltaOM=DeltaOm*2*pi*1e-6;
    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end
end

% Consistency enforcement
function grumble(parameters,DeltaOm)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isnumeric(DeltaOm)
    error('DeltaOm must be numeric.');
end
end
