% calculate the propagator of the CP pulse with Bterm
%
% input parameters:
% parameters: structure containing simulation parameters
% expt: the Map containing the experimental parameters
% v_CP: CP frequency
% H_IN: acting Hamiltonian without RF term
% Ix: the Ix spin Operator dictionary
% t: pulse length
% Ni_ENDOR: number of ENDOR nuclei
% N_spinSys: number of spin systems
% spin_system: Spinach spin system object
%
% output parameters:
% U: propagator for the pulse
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function U=kehl_cp_bterm(parameters,expt,v_CP,H_IN,Ix,t,Ni_ENDOR,N_spinSys,spin_system)
    if nargin<9
        spin_system=[];
    end


    % Check consistency
    grumble(parameters,expt,v_CP,H_IN,Ix,t,Ni_ENDOR,N_spinSys,spin_system);
    t_stepCP=1/(v_CP*parameters.N_stepRF);
    U_SL=eye(size(H_IN));
    H_SL=H_IN;

    for ll=1:parameters.N_stepRF
        % build the Hamiltonian for the increment
        if N_spinSys==1
            for mm=1:Ni_ENDOR
                H_SL=H_SL+2*expt("CP")*Ix{mm}*cos(2*pi*v_CP*t_stepCP*(ll-1));
            end
        else
            mm=1;
            H_SL=H_SL+2*expt("CP")*Ix{mm}*cos(2*pi*v_RF_CP*t_stepCP*(ll-1));
        end
        % create propagator and multiply with previous one
        if isempty(spin_system)
            U_step=expm(-1i*H_SL*t_stepCP);
        else
            U_step=full(propagator(spin_system,sparse(H_SL),t_stepCP));
        end
        U_SL=U_step*U_SL;
    end
    % build full propagator
    U=(U_SL^(t*v_CP));
end

function grumble(parameters,expt,v_CP,H_IN,Ix,t,Ni_ENDOR,N_spinSys,spin_system)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isnumeric(v_CP)
    error('v_CP must be numeric.');
end
if ~isnumeric(H_IN)
    error('H_IN must be numeric.');
end
if (~iscell(Ix))&&(~isnumeric(Ix))
    error('Ix must be numeric, or a cell array.');
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

