% Cross-polarisation frequency offset for Kehl ENDOR kernels. Syntax:
%
%      [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR_IN,...
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR_IN   - map containing ENDOR parameters.
%   v_off_S          - electron offset frequency.
%   HF_zz,NQI_zz     - effective hyperfine and quadrupole couplings.
%   m,kk             - nucleus and offset indices.
%
% Outputs:
%
%   v_CP             - cross-polarisation frequency.
%   paramsENDOR      - updated ENDOR parameter map.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_cp_offset.m>

function [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR_IN,...
        v_off_S,HF_zz,NQI_zz,m,kk)

    % Check consistency
    grumble(parameters,paramsENDOR_IN,v_off_S,HF_zz,NQI_zz,m,kk);

    % Get the explicit nuclear spin quantum numbers
    paramsENDOR=paramsENDOR_IN;
    spin_numbers=parameters.endor_spin_numbers;

    % Calculate the CP offset without the RF B term
    if parameters.Bterm==false
        if spin_numbers(m)==1/2
            if kk<0
                nu_alpha=sqrt((v_off_S+HF_zz(m)/2)^2+...
                    (parameters.electron_nutation/(2*pi))^2);
                nu_beta=sqrt((v_off_S-HF_zz(m)/2)^2+...
                    (parameters.electron_nutation/(2*pi))^2);

                if HF_zz(m)>0
                    offset_CP=1/2*(nu_alpha+nu_beta);
                else
                    offset_CP=1/2*(nu_alpha-nu_beta);
                end
            else
                nu_alpha=sqrt((v_off_S+HF_zz(m)/2)^2+...
                    (parameters.electron_nutation/(2*pi))^2);
                nu_beta=sqrt((v_off_S-HF_zz(m)/2)^2+...
                    (parameters.electron_nutation/(2*pi))^2);

                if HF_zz(m)>0
                    offset_CP=1/2*(nu_alpha-nu_beta);
                else
                    offset_CP=-1/2*(nu_alpha+nu_beta);
                end
            end

            % Select the opposite CP branch when requested
            if parameters.sel_CP==2
                offset_CP=-offset_CP;
            end
        elseif spin_numbers(m)==1

            % Legacy branch for spin-1 nuclei
            if kk==1
                v_off_S=freq_EPR(2)-freq_EPR(3);
            elseif kk==0
                v_off_S=0;
            elseif kk==-1
                v_off_S=freq_EPR(2)-freq_EPR(1);
            end

            omega_alpha=sqrt((v_off_S+HF_zz(m))^2+...
                (parameters.electron_nutation/(2*pi))^2);
            omega_beta=sqrt((v_off_S)^2+...
                (parameters.electron_nutation/(2*pi))^2);
            omega_gamma=sqrt((v_off_S-HF_zz(m))^2+...
                (parameters.electron_nutation/(2*pi))^2);

            if kk==1
                if parameters.sel_CP==3
                    offset_CP=-(-omega_beta+omega_alpha+3*NQI_zz(m))/2;
                elseif parameters.sel_CP==1
                    offset_CP=-(omega_beta-omega_alpha+3*NQI_zz(m))/2;
                end
            elseif kk==0
                if parameters.sel_CP==4
                    offset_CP=(omega_gamma-omega_beta+3*NQI_zz(m))/2;
                elseif parameters.sel_CP==3
                    offset_CP=-(-omega_beta-omega_alpha+3*NQI_zz(m))/2;
                elseif parameters.sel_CP==2
                    offset_CP=(-omega_gamma+omega_beta+3*NQI_zz(m))/2;
                elseif parameters.sel_CP==1
                    offset_CP=-(+omega_beta+omega_alpha+3*NQI_zz(m))/2;
                end
            elseif kk==(-1)
                if parameters.sel_CP==2
                    offset_CP=(omega_gamma+omega_beta+3*NQI_zz(m))/2;
                elseif parameters.sel_CP==1
                    offset_CP=(-omega_gamma-omega_beta+3*NQI_zz(m))/2;
                end
            end
        end
    else

        % Use externally supplied CP start when RF B term is active
        offset_CP=parameters.cp_start_hz;
    end

    % Write the selected CP frequency into the ENDOR parameter map
    v_L=paramsENDOR('v_L');
    v_CP=v_L(m)-offset_CP;
    paramsENDOR('yCoords')=v_CP;

end

% Consistency enforcement
function grumble(parameters,paramsENDOR_IN,v_off_S,HF_zz,NQI_zz,m,kk)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isa(paramsENDOR_IN,'containers.Map')
        error('paramsENDOR_IN must be a containers.Map object.');
    end
    if ~isnumeric(v_off_S)
        error('v_off_S must be numeric.');
    end
    if ~isnumeric(HF_zz)
        error('HF_zz must be numeric.');
    end
    if ~isnumeric(NQI_zz)
        error('NQI_zz must be numeric.');
    end
    if ~isnumeric(m)
        error('m must be numeric.');
    end
    if ~isnumeric(kk)
        error('kk must be numeric.');
    end
end

