%KEHL_LINE_BROADEN Line broadening for Kehl ENDOR spectra.
%
%   Spinach architecture migration May 2026 Talos

function data_conv=kehl_line_broaden(data,parameters)

    % Check consistency
    grumble(data,parameters);

    % Get ENDOR sweep width
    sw=parameters.endor_range_hz;
    if parameters.Lorentzian==1

        % Lorentian
        Deltaend_L=parameters.lw_L*pi/(sw);
        endintens1=ifft(data(:));
        for i=1:size(data,2)

            % Lorentzian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_L*i);
        end
        endintens(1)=0.37*endintens(1);
        endamp_2d_L_conv(:)=real(fft(endintens));

        if parameters.Gaussian==1
            % Gaussian
            lw=parameters.lw_G;
            Deltaend_G=(lw*pi/(sw*sqrt(2*log(2))))^2;
            endintens1=ifft(endamp_2d_L_conv(:)+1000);
            for i=1:size(data,2)

                % Gaussian line shape
                endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
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

        %size(data,2)
        for i=1:length(data)

            % Gaussian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
        end
        endintens(1)=0.5*endintens(1);
        data_conv(:)=-abs(fft(endintens(:)));

    else
        data_conv=data;
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
