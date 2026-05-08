% start CP step calculation
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% expt: the Map containing the experimental parameters
% opt: the Map containing the optional paramters
%
% output parameters:
% rho: the final spin density matrix
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function [rho]=kehl_cp_step(constants,spinSys,spinOps,expt,opt)

    % Check consistency
    grumble(constants,spinSys,spinOps,expt,opt);
    if opt("Relax")==0
        opt("Relax_step")=0;
        opt("Relax")=1;
        opt("temp_eff")=1;
    elseif ~isKey(opt,"Relax_step")
        opt("Relax_step")=1;
    end

    % prepare EPR and ENDOR calculations
    paramsEPR=kehl_prep_epr(spinSys,expt);
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);


    % for SC no 2D scan, and for Powder define y-axis
    if opt("powder")==false
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
    if opt('freqDomain')==false
        epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
    else
        epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
    end

    % ENDOR calculation

    % starts the actual calculation
    rho=kehl_cp_step_calc(constants,spinOps,spinSys,expt,opt,paramsENDOR,epr);

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

