% start CP ENDOR calculation
% is directly called from the set up
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
%
% output parameters:
% endor_amp: endor amplitude values
% endor_amp_conv: endor values convoluted with lb
% x_coords: coords of rf axis
% y_coords: coords of cp axis
% v_L: nuclear Larmor frequency
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% prepare EPR and ENDOR calculations

function [endor_amp,endor_amp_conv,x_coords,y_coords,v_L]=kehl_cp_sim(constants,spinSys,spinOps,expt,parameters)

    % Check consistency
    grumble(constants,spinSys,spinOps,expt,parameters);
    parameters=kehl_defaults(parameters);
    paramsEPR=kehl_prep_epr(spinSys,expt);
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);

    x_coords=paramsENDOR("x_coords");
    x_coords=x_coords(:,1)';
    v_L=paramsENDOR("v_L");

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
            % step_CP=0;
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
    if parameters.Relax==true

        % starts the actual calculation with relaxation
        endor_amp=kehl_cp_calc_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,epr);
    else

        % starts the actual calculation
        endor_amp=kehl_cp_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,epr);
    end

    if parameters.powder==false
        y_coords=1;
    end

    % convolution with lb
    endor_amp_conv=kehl_line_broaden(endor_amp,parameters,expt("range_EN"));

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

