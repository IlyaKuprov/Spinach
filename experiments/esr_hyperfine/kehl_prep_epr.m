% calculates/sets necessary parameters for the EPR calculation from
% the experimental values
%
% input parameters:
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
%
% output parameters:
% paramsEPR: the Map containing the EPR parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)



function paramsEPR=kehl_prep_epr(spinSys,expt)

    % Check consistency
    grumble(spinSys,expt);
    g_iso=spinSys("g_iso");
    obsField=expt("Field");

    gBepr=obsField*g_iso;
    geff=g_iso;
    B=gBepr/geff;
    veff=geff*obsField*9.27401*1e-24/(6.62607*1e-34);

    % EPR simulation range
    fieldCenter=gBepr/g_iso;
    fieldmin=fieldCenter-0.260;
    fieldmax=fieldCenter+0.260;

    % EPR x-axis definition

    % no of points in Field dimension
    Npts=(round((fieldmax-fieldmin)/expt("deltaField"))+1);
    field=linspace(fieldmin,fieldmax,Npts);

    paramsEPR=containers.Map;
    paramsEPR("Npts")=Npts;
    paramsEPR("fieldAxis")=field;

end

function grumble(spinSys,expt)
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end

