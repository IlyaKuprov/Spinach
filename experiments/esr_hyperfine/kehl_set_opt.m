% sets the option parameters in the Map opt to be used in the
% simulation routine
%
% input parameters:
% spinSys: the Map describing the spin system
% defaults: the Map with the new values for the options
%
% output parameters:
% opt: Map with the new option values
%
% Febuary 2024 A. Kehl (akehl@gwdg.de)
%

function opt=kehl_set_opt(spinSys,values)

    % Check consistency
    grumble(spinSys,values);
    opt=containers.Map;

    % Is the powder pattern used?
    % default false
    if isKey(values,'powder')
        opt("powder")=values("powder");
    else
        opt("powder")=false;
    end

    % Is the EPR spectrum calculated in the freq. or field domain?
    % default field domain
    if isKey(values,'freqDomain')
        opt('freqDomain')=values('freqDomain');
    else
        opt('freqDomain')=false;
    end

    % Is the Bterm included?
    % default false
    if isKey(values,'use_Bterm')
        opt("Bterm")=values("use_Bterm");
    else
        opt("Bterm")=false;
    end

    % How many steps are calculated for the steps in the RF pulse when the
    % Bterm is included?
    % default 100
    if isKey(values,'N_stepRF')
        opt("N_stepRF")=values("N_stepRF");
    else
        opt("N_stepRF")=100;
    end

    % Is relaxation included?
    % default false
    if isKey(values,'Relax')
        opt('Relax')=values('Relax');
    else
        opt('Relax')=false;
    end

    % Are temperature effects used?
    % default false
    if isKey(values,'temp_eff')
        opt("temp_eff")=values("temp_eff");
    else
        opt("temp_eff")=false;
    end

    % Use temperature effects if relaxation is used.
    if isKey(values,"Relax")
        if values("Relax")==true
            if opt("temp_eff")==false
                warning('The temperature needs to be considered if relaxation is used.');
            end
            opt("temp_eff")=true;
        end
    end

    % Set the temperature.
    % default 10K
    if isKey(values,'T')
        opt("T")=values("T");
    else
        opt("T")=10;
    end

    % Choose the EPR line for single crystal calculation.
    % default: EPR alpha
    if isKey(values,'sel_I')
        opt("sel_I")=values("sel_I");
    else
        opt("sel_I")=1;
    end

    % Choose the CP condition for CP ENDOR.
    % default reduce intensity of line with lowest RF freq
    if isKey(values,'sel_CP')
        opt("sel_CP")=values("sel_CP");
    else
        opt("sel_CP")=1;
    end

    % How fine is the grid for the powder pattern?
    % default 20 (low resolution, but time efficient!)
    if isKey(values,'Nang')
        opt("Nang")=values("Nang");
    else
        opt("Nang")=20;
    end

    % How broad is the excitation bandwidth?
    % default 20kHz
    if isKey(values,'nwidth')
        opt("nwidth")=values("nwidth");
    else
        opt("nwidth")=20;
    end

    % Is Lorentian lb used?
    % default false
    if isKey(values,'Lorentzian')
        opt("Lorentzian")=values("Lorentzian");
    else
        opt("Lorentzian")=false;
    end

    % Lorentian lb value
    % default 20kHz
    if isKey(values,'lw_L')
        opt("lw_L")=values("lw_L")*1000;
    else
        opt("lw_L")=20*1e3;
    end

    % Is Gaussian lb used?
    % default false
    if isKey(values,'Gaussian')
        opt("Gaussian")=values("Gaussian");
    else
        opt("Gaussian")=false;
    end

    % Gaussian lb value
    % default 20kHz
    if isKey(values,'lw_G')
        opt("lw_G")=values("lw_G")*1000;
    else
        opt("lw_G")=20*1e3;
    end

end

function grumble(spinSys,values)
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(values,'containers.Map')
    error('values must be a containers.Map object.');
end
end

