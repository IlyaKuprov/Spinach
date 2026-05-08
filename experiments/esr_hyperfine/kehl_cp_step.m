% start CP step calculation
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
%
% output parameters:
% rho: the final spin density matrix
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function [rho]=kehl_cp_step(constants,spinSys,spinOps,expt,parameters)

    % Check consistency
    grumble(constants,spinSys,spinOps,expt,parameters);
    parameters=kehl_defaults(parameters);
    if parameters.Relax==0
        parameters.Relax_step=0;
        parameters.Relax=1;
        parameters.temp_eff=1;
    elseif ~isfield(parameters,"Relax_step")
        parameters.Relax_step=1;
    end

    % prepare EPR and ENDOR calculations
    paramsEPR=kehl_prep_epr(spinSys,expt);
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);


    % for SC no 2D scan, and for Powder define y-axis
    if parameters.powder==false
        expt("Npts_CP")=1;
        expt("range_CP")=0;
    else
        if expt("Npts_CP")>1
            step_CP=expt("range_CP")/(expt("Npts_CP")-1);
            y_coords=zeros(expt("Npts_CP"));
            for ii=1:expt("Npts_CP")
                y_coords(ii)=expt("start_CP")+(ii-1)*step_CP;
            end
            y_coords=y_coords(:,1)';
        else
            y_coords=expt("start_CP");
            y_coords=y_coords(:,1)';
        end
        paramsENDOR("y_coords")=y_coords;
    end

    % EPR calculation
    if parameters.freqDomain==false
        epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    else
        epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    end

    % ENDOR calculation

    % starts the actual calculation
    rho=kehl_cp_step_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,epr);

end

function grumble(constants,spinSys,spinOps,expt,parameters)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

