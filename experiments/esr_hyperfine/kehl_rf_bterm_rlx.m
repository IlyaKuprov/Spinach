%KEHL_RF_BTERM_RLX Relaxation-space RF B-term propagator.
%
%   U=KEHL_RF_BTERM_RLX(PARAMETERS,V_RF,HFREE,IY,T,N_ENDOR,
%   N_SPIN_SYSTEMS,R,SPIN_SYSTEM) creates the propagator for an RF pulse
%   with the RF B term when relaxation is included.
%
%   Inputs:
%
%      PARAMETERS     - Kehl ENDOR context parameters.
%      V_RF           - RF frequency.
%      HFREE          - Hamiltonian without the RF term.
%      IY             - cell array of nuclear Iy spin operators.
%      T              - pulse length.
%      N_ENDOR        - number of ENDOR nuclei.
%      N_SPIN_SYSTEMS - number of spin systems.
%      R              - relaxation superoperator.
%      SPIN_SYSTEM    - Spinach spin system structure.
%
%   Output:
%
%      U - propagator for the pulse.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function U=kehl_rf_bterm_rlx(parameters,v_RF,Hfree,Iy,t,n_endor,...
                             n_spin_systems,R,spin_system)
if nargin<9
    spin_system=[];
end

% Check consistency
grumble(parameters,v_RF,Hfree,Iy,t,n_endor,n_spin_systems,R,spin_system);

% Incrementation calculation for the RF pulse
t_stepRF=1/(v_RF*parameters.N_stepRF);
U_RF=eye(size(R));

for n=1:parameters.N_stepRF
    if n_spin_systems==1
        HRF=Hfree;
        for k=1:n_endor
            HRF=HRF-2*parameters.nuclear_nutation*Iy{k}*...
                cos(2*pi*v_RF*t_stepRF*(n-1));
        end
    else
        HRF=Hfree-2*parameters.nuclear_nutation*Iy{1}*...
            cos(2*pi*v_RF*t_stepRF*(n-1));
    end
    G=R-1i*full(hilb2liouv(sparse(HRF),'comm'));
    if isempty(spin_system)
        U_step=expm(G*t_stepRF);
    else
        U_step=full(propagator(spin_system,1i*sparse(G),t_stepRF));
    end
    U_RF=U_step*U_RF;
end

U=U_RF^(t*v_RF);

end

% Consistency enforcement
function grumble(parameters,v_RF,Hfree,Iy,t,n_endor,n_spin_systems,R,...
                 spin_system)
if ~isstruct(parameters)
    error('parameters must be a structure.');
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
if ~isnumeric(R)
    error('R must be numeric.');
end
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
   (~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
end
