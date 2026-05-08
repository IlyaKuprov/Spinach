% creates the propagator for the RF pulse if the Bterm is considered
%
% input parameters:
% parameters: structure containing simulation parameters
% expt: the Map containing the experimental parameters
% v_RF: RF frequency
% Hfree: acting Hamiltonian without RF term
% Iy: cell array of nuclear Iy spin operators
% t: pulse length
% n_endor: number of ENDOR nuclei
% n_spin_systems: number of spin systems
% spin_system: Spinach spin system object
%
% output parameters:
% U: propagator for the pulse
%
% February 2024 A. Kehl (akehl@gwdg.de)


function U=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t,n_endor,n_spin_systems,spin_system)
    if nargin<9
        spin_system=[];
    end


    % Check consistency
    grumble(parameters,expt,v_RF,Hfree,Iy,t,n_endor,n_spin_systems,spin_system);
    % Incrementation calculation for the RF pulse
    t_stepRF=1/(v_RF*parameters.N_stepRF);
    U_RF=eye(size(Hfree));

    for ll=1:parameters.N_stepRF
        if n_spin_systems==1
            HRF=Hfree;
            for mm=1:n_endor
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

function grumble(parameters,expt,v_RF,Hfree,Iy,t,n_endor,n_spin_systems,spin_system)
if ~isstruct(parameters)
    error('parameters must be a structure.');
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
if ~isnumeric(n_endor)
    error('n_endor must be numeric.');
end
if ~isnumeric(n_spin_systems)
    error('n_spin_systems must be numeric.');
end
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
end

