% calculate the cp offset to hit one cp transition for the sc
% calculation
% only for I=1/2 and I=1
%
% input parameters:
% parameters: structure containing simulation parameters
% paramsENDOR_IN: the Map containing the ENDOR parameters
% expt: the Map containing the experimental parameters
% spinSys: the Map describing the spin system
% v_off_S: electron offset freq
% HF_zz: effective HF coupling value
% NQI_zz: effective NQC value
% m: number of the nucleus (loop iteration)
% kk: number of the offset (loop iteration)
%
% output parameters:
% v_CP: cp frequency
% paramsENDOR: updated Map
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR_IN,expt,spinSys,v_off_S,HF_zz,NQI_zz,m,kk)

    % Check consistency
    grumble(parameters,paramsENDOR_IN,expt,spinSys,v_off_S,HF_zz,NQI_zz,m,kk);
    paramsENDOR=paramsENDOR_IN;
    I=spinSys("I");

    if parameters.Bterm==false
         if I(m)==1/2
            if kk<0

                nu_alpha=sqrt((v_off_S+HF_zz(m)/2)^2+(expt("oneE")/(2*pi))^2);
                nu_beta=sqrt((v_off_S-HF_zz(m)/2)^2+(expt("oneE")/(2*pi))^2);

                if HF_zz(m)>0
                    offset_CP=1/2*(nu_alpha+nu_beta);
                else
                    offset_CP=1/2*(nu_alpha-nu_beta);
                end

            else

                nu_alpha=sqrt((v_off_S+HF_zz(m)/2)^2+(expt("oneE")/(2*pi))^2);
                nu_beta=sqrt((v_off_S-HF_zz(m)/2)^2+(expt("oneE")/(2*pi))^2);

                if HF_zz(m)>0
                    offset_CP=1/2*(nu_alpha-nu_beta);
                else
                    offset_CP=-1/2*(nu_alpha+nu_beta);
                end
            end

            if parameters.sel_CP==2
                offset_CP=-offset_CP;
            end


        elseif I(m)==1
          % from diagonalization for 2D
            if kk==1
                v_off_S=freq_EPR(2)-freq_EPR(3);
            elseif kk==0
                v_off_S=0;
            elseif kk==-1
                v_off_S=freq_EPR(2)-freq_EPR(1);
            end


            omega_alpha=sqrt((v_off_S+HF_zz(m))^2+(expt("oneE")/(2*pi))^2);
            omega_beta=sqrt((v_off_S)^2+(expt("oneE")/(2*pi))^2);
            omega_gamma=sqrt((v_off_S-HF_zz(m))^2+(expt("oneE")/(2*pi))^2);


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
        offset_CP=expt("start_CP");
    end
    v_L=paramsENDOR("v_L");
    v_CP=v_L(m)-offset_CP;
    paramsENDOR("yCoords")=v_CP;
end

function grumble(parameters,paramsENDOR_IN,expt,spinSys,v_off_S,HF_zz,NQI_zz,m,kk)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR_IN,'containers.Map')
    error('paramsENDOR_IN must be a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
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

