% Microwave pulse excitation factor for field-domain selection. Syntax:
%
%      scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   DeltaB           - field offset.
%
% Outputs:
%
%   scalefactor      - orientation excitation factor.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ori_field_pulsescale.m>

function scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB)

    % Check consistency
    grumble(parameters,DeltaB);

    % Load microwave pulse shape
    pulse=fopen(parameters.pulse_file);
    if pulse<0
        error('parameters.pulse_file could not be opened.');
    end
    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose(pulse);

    % Compute normalised excitation profile
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));
    n=1000;
    pulse_y=pulse_data{3}+1i*pulse_data{4};
    Y=abs(fftshift(fft(pulse_y,(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);

    % Select multipulse excitation profile if requested
    if isfield(parameters,'multipulses')&&parameters.multipulses==true
        Y=y;
    end

    % Convert field offset to pulse-profile bins
    DeltaOM=DeltaB*parameters.constants('GE')*1e-6;

    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0&&binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end

end

% Consistency enforcement
function grumble(parameters,DeltaB)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if (~isfield(parameters,'pulse_file'))||(~ischar(parameters.pulse_file)&&~isstring(parameters.pulse_file))
        error('parameters.pulse_file must be a character string.');
    end
    if ~isnumeric(DeltaB)
        error('DeltaB must be numeric.');
    end
    if (~isfield(parameters,'constants'))||(~isa(parameters.constants,'containers.Map'))
        error('parameters.constants must be a containers.Map object.');
    end
end

