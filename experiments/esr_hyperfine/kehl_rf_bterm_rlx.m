% Relaxation-space RF B-term propagator for Kehl ENDOR kernels. Syntax:
%
%      U=kehl_rf_bterm_rlx(parameters,v_RF,Hfree,Iy,t,n_endor,R,...
%                          spin_system)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   v_RF             - radiofrequency offset.
%   Hfree            - Hamiltonian without the RF term.
%   Iy               - cell array of nuclear Iy spin operators.
%   t                - pulse length.
%   n_endor          - number of ENDOR nuclei.
%   R                - relaxation superoperator.
%   spin_system      - Spinach spin system structure.
%
% Outputs:
%
%   U                - RF pulse propagator.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_rf_bterm_rlx.m>

function U=kehl_rf_bterm_rlx(parameters,v_RF,Hfree,Iy,t,n_endor,R,...
        spin_system)

    % Check consistency
    grumble(parameters,v_RF,Hfree,Iy,t,n_endor,R,spin_system);

    % Incrementation calculation for the RF pulse
    t_stepRF=1/(v_RF*parameters.N_stepRF);
    U_RF=eye(size(R));

    for n=1:parameters.N_stepRF
        HRF=Hfree;
        for k=1:n_endor
            HRF=HRF-2*parameters.nuclear_nutation*Iy{k}*...
                cos(2*pi*v_RF*t_stepRF*(n-1));
        end
        G=R-1i*full(hilb2liouv(sparse(HRF),'comm'));
        U_step=full(propagator(spin_system,1i*sparse(G),t_stepRF));
        U_RF=U_step*U_RF;
    end

    U=U_RF^(t*v_RF);

end

% Consistency enforcement
function grumble(parameters,v_RF,Hfree,Iy,t,n_endor,R,spin_system)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isnumeric(v_RF)
        error('v_RF must be numeric.');
    end
    if ~isnumeric(Hfree)
        error('Hfree must be numeric.');
    end
    if ~iscell(Iy)
        error('Iy must be a cell array.');
    end
    if ~isnumeric(t)
        error('t must be numeric.');
    end
    if ~isnumeric(n_endor)
        error('n_endor must be numeric.');
    end
    if ~isnumeric(R)
        error('R must be numeric.');
    end
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
end

