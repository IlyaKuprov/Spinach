% convolutes the input data with a Lorentian and/or a
% Gaussian line broadening
% CAUTION: this is a beta version and not fully tested
% convolution for time domain traces can easily be performed in Origin,
% etc. by multiplication with the corresponding exponential function
%
% input parameters:
% data: data to be convoluted
% opt: options Map, should contain keys ('Lorentian'/'Gaussian' and
%       'lw_L'/'lw_G')
%
% output parameters:
% data_conv: convoluted data
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function data_conv=kehl_time_lb(data,opt)

    % Check consistency
    grumble(data,opt);
data_n=(data-min(data))/(max(data)-min(data))-0.5;

if opt("Lorentzian")==1
        warning('This is a beta version and has not been fully tested.');


        Deltaend_L=opt('lw_L')*pi/0.1;
        % Lorentian
        endintens1=(data_n(:));
        for i=1:size(data,2)

            % Lorentzian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_L*i);
        end
        endintens(1)=0.37*endintens(1);
        endamp_2d_L_conv(:)=real((endintens));

        if opt("Gaussian")==1
            lw=opt('lw_G');
            Deltaend_G=(lw*pi/(0.1*sqrt(2*log(2))))^2;
            % Gaussian
            endintens1=(endamp_2d_L_conv(:));
            for i=1:size(data_n,2)

                % Gaussian line shape
                endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
            end
            endintens(1)=0.5*endintens(1);
            data_conv(:)=-abs((endintens));

        else
            data_conv=endamp_2d_L_conv(:);
        end
    elseif opt("Gaussian")==1
        warning('This is a beta version and has not been fully tested.');

        lw=opt('lw_G');
        Deltaend_G=(lw*pi/(0.1*sqrt(2*log(2))))^2;

        % Gaussian
        data=data+1000;
        endintens1=(data_n(:));
        for i=1:size(data_n,2)

            % Gaussian line shape
            endintens(i)=endintens1(i)*exp(-Deltaend_G*i^2/2);
        end
        endintens(1)=0.5*endintens(1);
        data_conv(:)=-abs((endintens));

    else
        data_conv=data;
end
end

function grumble(data,opt)
if (~isnumeric(data))&&(~ischar(data))&&(~isstring(data))
    error('data must be numeric, a character string, or a string scalar.');
end
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
end

