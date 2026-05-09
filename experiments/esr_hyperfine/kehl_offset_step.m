% Offset integration step for Kehl ENDOR kernels. Syntax:
%
%      timestep=kehl_offset_step(offset,reference_offset,n_steps)
%
% Parameters:
%
%   offset           - current electron angular-frequency offset.
%   reference_offset - first electron angular-frequency offset.
%   n_steps          - integration divisor.
%
% Outputs:
%
%   timestep         - finite integration timestep.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_offset_step.m>

function timestep=kehl_offset_step(offset,reference_offset,n_steps)

    % Check consistency
    grumble(offset,reference_offset,n_steps);

    % Use the current offset when it is non-zero
    if offset~=0
        timestep=2*pi/(offset*n_steps);
        return
    end

    % Fall back to the reference offset if available
    if reference_offset~=0
        timestep=2*pi/(reference_offset*n_steps);
        return
    end

    % Degenerate offset manifolds do not need oscillation integration
    timestep=0;

end

% Consistency enforcement
function grumble(offset,reference_offset,n_steps)
    if (~isnumeric(offset))||(~isscalar(offset))||(~isfinite(offset))
        error('offset must be a finite numeric scalar.');
    end
    if (~isnumeric(reference_offset))||(~isscalar(reference_offset))||...
            (~isfinite(reference_offset))
        error('reference_offset must be a finite numeric scalar.');
    end
    if (~isnumeric(n_steps))||(~isscalar(n_steps))||(~isfinite(n_steps))||...
            (n_steps<=0)
        error('n_steps must be a positive finite numeric scalar.');
    end
end

