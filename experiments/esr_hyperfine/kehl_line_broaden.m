% Line broadening for Kehl ENDOR spectra. Syntax:
%
%      data_conv=kehl_line_broaden(data,parameters)
%
% Parameters:
%
%   data             - ENDOR amplitude array.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   data_conv        - line-broadened ENDOR amplitude array.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_line_broaden.m>

function data_conv=kehl_line_broaden(data,parameters)

    % Check consistency
    grumble(data,parameters);

    % Return immediately when broadening is disabled
    if (parameters.Lorentzian~=1)&&(parameters.Gaussian~=1)
        data_conv=data;
        return
    end

    % Get ENDOR sweep width
    sw=parameters.endor_range;
    if parameters.Lorentzian==1

        % Lorentzian broadening
        Deltaend_L=parameters.lw_L*pi/(sw);
        endintens1=ifft(data(:));
        for point_idx=1:size(data,2)

            % Lorentzian line shape
            endintens(point_idx)=endintens1(point_idx)*exp(-Deltaend_L*point_idx);
        end
        endintens(1)=0.37*endintens(1);
        endamp_2d_L_conv(:)=real(fft(endintens));

        if parameters.Gaussian==1
            % Gaussian
            lw=parameters.lw_G;
            Deltaend_G=(lw*pi/(sw*sqrt(2*log(2))))^2;
            endintens1=ifft(endamp_2d_L_conv(:)+1000);
            for point_idx=1:size(data,2)

                % Gaussian line shape
                endintens(point_idx)=endintens1(point_idx)*exp(-Deltaend_G*point_idx^2/2);
            end
            endintens(1)=0.5*endintens(1);
            data_conv(:)=-abs(fft(endintens));

        else
            data_conv=endamp_2d_L_conv(:);
        end
    elseif parameters.Gaussian==1
        % Gaussian
        lw=parameters.lw_G;
        Deltaend_G=(lw*pi/(sw*sqrt(2*log(2))))^2;
        data=data+1000;
        endintens1=ifft(data(:));
        for point_idx=1:length(data)

            % Gaussian line shape
            endintens(point_idx)=endintens1(point_idx)*exp(-Deltaend_G*point_idx^2/2);
        end
        endintens(1)=0.5*endintens(1);
        data_conv(:)=-abs(fft(endintens(:)));

    end
end
function grumble(data,parameters)
    if ~isnumeric(data)
        error('data must be numeric.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

