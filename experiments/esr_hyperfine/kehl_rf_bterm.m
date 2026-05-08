% creates the propagator for the RF pulse if the Bterm is considered
%
% input parameters:
% opt: the Map containing the optional paramters
% expt: the Map containing the experimental parameters
% v_RF: RF frequency
% Hfree: acting Hamiltonian without RF term
% Iy: the Iy spin Operator dictionary
% t: pulse length
% Ni_ENDOR: number of ENDOR nuclei
% N_spinSys: number of spin systems
% spin_system: Spinach spin system object
%
% output parameters:
% U: propagator for the pulse
%
% February 2024 A. Kehl (akehl@gwdg.de)


function U=kehl_rf_bterm(opt,expt,v_RF,Hfree,Iy,t,Ni_ENDOR,N_spinSys,spin_system)
    if nargin<9
        spin_system=[];
    end


    % Check consistency
    grumble(opt,expt,v_RF,Hfree,Iy,t,Ni_ENDOR,N_spinSys,spin_system);
    % Incrementation calculation for the RF pulse
    t_stepRF=1/(v_RF*opt("N_stepRF"));
    U_RF=eye(size(Hfree));

    for ll=1:opt("N_stepRF")
        if N_spinSys==1
            HRF=Hfree;
            for mm=1:Ni_ENDOR
                HRF=HRF+2*(expt("oneN"))*Iy{mm}*cos(2*pi*v_RF*t_stepRF*(ll-1));
            end
        else
            HRF=Hfree+2*(expt("oneN"))*Iy{1}*cos(2*pi*v_RF*t_stepRF*(ll-1));
        end
        if isempty(spin_system)
            U_step=expm(-1i*HRF*t_stepRF);
        else
            U_step=full(propagator(spin_system,sparse(HRF),t_stepRF));
        end
        U_RF=U_step*U_RF;
    end

    U=(U_RF^(t*v_RF));

end

function grumble(opt,expt,v_RF,Hfree,Iy,t,Ni_ENDOR,N_spinSys,spin_system)
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isnumeric(v_RF)
    error('v_RF must be numeric.');
end
if ~isnumeric(Hfree)
    error('Hfree must be numeric.');
end
if (~iscell(Iy))&&(~isnumeric(Iy))
    error('Iy must be numeric, or a cell array.');
end
if ~isnumeric(t)
    error('t must be numeric.');
end
if ~isnumeric(Ni_ENDOR)
    error('Ni_ENDOR must be numeric.');
end
if ~isnumeric(N_spinSys)
    error('N_spinSys must be numeric.');
end
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
end

