%KEHL_ORI_FIELD_PULSESCALE Microwave pulse excitation factor for field-domain selection.
%
%   Spinach architecture migration May 2026 Talos

function scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB,CONST1)

% Check consistency
grumble(parameters,DeltaB,CONST1);

    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % parameters: Kehl ENDOR context parameters
    % DeltaB: offset of the actual field from the effective resonance field
    % CONST1: constant to convert between field and frequency
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    data=parameters.pulse_file;
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y,(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isfield(parameters,'multipulses')&&parameters.multipulses==true
       Y=y;
    end

    %Delta omega/2pi in MHz
    DeltaOM=DeltaB*CONST1*10^4*2*pi;

    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end

end

% Consistency enforcement
function grumble(parameters,DeltaB,CONST1)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isnumeric(DeltaB)
    error('DeltaB must be numeric.');
end
if ~isnumeric(CONST1)
    error('CONST1 must be numeric.');
end
end
