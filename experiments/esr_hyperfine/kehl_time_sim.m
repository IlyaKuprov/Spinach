% start time domain  ENDOR calculation
% is directly called from the set up
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% expt: the Map containing the experimental parameters
% opt: the Map containing the optional paramters
%
% output parameters:
% endor_amp: endor amplitude values
% endor_amp_conv: endor values convoluted with lb
% x_coords: coords of rf axis
% v_L: nuclear Larmor frequency
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function [endor_amp,endor_amp_conv,x_coords,v_L]=kehl_time_sim(constants,spinSys,spinOps,expt,opt)

    % Check consistency
    grumble(constants,spinSys,spinOps,expt,opt);
    t=expt('t');
    expt('t_start')=t(6);

    % prepare EPR and ENDOR calculations
    paramsEPR=kehl_prep_epr(spinSys,expt);
    paramsENDOR=kehl_prep_time(constants,spinSys,expt);

    x_coords=paramsENDOR("x_coords");
    x_coords=x_coords(:,1)';
    v_L=paramsENDOR("v_L");


    % EPR calculation
    if opt('freqDomain')==false
        epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
    else
        epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
    end

    % ENDOR calculation
    if opt("Relax")==true

        % starts the actual calculation with relaxation
        endor_amp=kehl_time_rlx(constants,spinOps,spinSys,expt,opt,paramsENDOR,epr);
    else

        % starts the actual calculation
        endor_amp=kehl_time_calc(constants,spinOps,spinSys,expt,opt,paramsENDOR,epr);
    end

    % convolution with lb
    endor_amp_conv=kehl_time_lb(endor_amp,opt);
end

function grumble(constants,spinSys,spinOps,expt,opt)
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
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
end

