% calculates/sets necessary parameters for the ENDOR calculation from
% the experimental values
%
% input parameters:
% constants:the Map containing the constants
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
%
% output parameters:
% paramsENDOR: the Map containing the ENDOR parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)

% Larmor frequencies

function paramsENDOR=kehl_prep_endor(constants,spinSys,expt)

    % Check consistency
    grumble(constants,spinSys,expt);
    Nuclei=spinSys("Nuclei");
    v_L=zeros(size(Nuclei,2),1);
    for i=1:size(Nuclei,2)
        v_L(i)=kehl_nuc_gamma(constants,Nuclei{i})*expt("Field");
    end

    % number of points and values for x-axis
    Npts_EN=(round(expt("range_EN")/expt("res_EN")));
    if isKey(expt,"RF_start")

        % ENDOR start x-axis
        start_EN=expt("RF_start");
    else

        % ENDOR start x-axis
        start_EN=v_L(1)-expt("range_EN")/2;
    end

    % X-Axis steps
    step_EN=expt("range_EN")/(Npts_EN-1);
    x_coords=zeros(Npts_EN);
    for ii=1:Npts_EN

        %-v_L;
        x_coords(ii)=start_EN+(ii-1)*step_EN;
    end

    paramsENDOR=containers.Map;

    paramsENDOR("v_L")=v_L;
    paramsENDOR("start_EN")=start_EN;
    paramsENDOR("step_EN")=step_EN;
    paramsENDOR("range_EN")=expt("range_EN");
    paramsENDOR("Npts_EN")=Npts_EN;
    paramsENDOR("x_coords")=x_coords;


end

function grumble(constants,spinSys,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end

