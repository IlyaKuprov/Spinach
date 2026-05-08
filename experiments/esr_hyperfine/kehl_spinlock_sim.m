% start Spinlock ENDOR calculation
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
% v_L: nuclear Larmor frequency
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% prepare EPR and ENDOR calculations

function [endor_amp,endor_amp_conv,x_coords,v_L]=kehl_spinlock_sim(constants,spinSys,spinOps,expt,parameters)

    % Check consistency
    grumble(constants,spinSys,spinOps,expt,parameters);
    parameters=kehl_defaults(parameters);
    paramsEPR=kehl_prep_epr(spinSys,expt);
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);

    x_coords=paramsENDOR("x_coords");
    x_coords=x_coords(:,1)';
    v_L=paramsENDOR("v_L");

    % EPR calculation
    if parameters.freqDomain==false
        epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    else
        epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    end

    % ENDOR calculation
    if parameters.Relax==true

        % starts the actual calculation with relaxation
        endor_amp=kehl_spinlock_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,epr);
    else

        % starts the actual calculation
        endor_amp=kehl_spinlock_calc(constants,spinOps,spinSys,expt,parameters,paramsENDOR,epr);
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

