% Electron spin offsets for Kehl ENDOR kernels. Syntax:
%
%      offsets=kehl_offsets(parameters,spin_system,B,euler_angles)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   spin_system      - Zeeman-Liouville Spinach spin system.
%   B                - magnetic field in Tesla.
%   euler_angles     - Kehl orientation angles stored by the context.
%
% Outputs:
%
%   offsets          - electron spin angular-frequency offsets for all selected EPR transitions.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_offsets.m>

function offsets=kehl_offsets(parameters,spin_system,B,euler_angles)

    % Check consistency
    grumble(parameters,spin_system,B,euler_angles);

    % Get Liouville-space EPR transition angular frequencies
    [freq_epr,moments]=kehl_epr_transitions(spin_system,parameters,...
        euler_angles,'frequency',B);

    % Keep microwave-active transitions
    if isempty(freq_epr)
        offsets=0;
        return
    end
    hit_list=moments>max(moments)*1e-10;
    freq_epr=freq_epr(hit_list);
    if isempty(freq_epr)
        offsets=0;
        return
    end

    % Sort transition angular frequencies for offset construction
    freq_epr=sort(freq_epr(:).');
    offsets=zeros(1,numel(freq_epr));

    % Calculate offsets around the transition centre
    if mod(numel(freq_epr),2)==0
        for n=1:numel(freq_epr)/2
            offsets(n)=-(freq_epr(numel(freq_epr)+1-n)-freq_epr(n))/2;
        end
        for n=numel(freq_epr)/2+1:numel(freq_epr)
            offsets(n)=-offsets(numel(freq_epr)+1-n);
        end
    else
        centre=numel(freq_epr)/2+0.5;
        for n=1:numel(freq_epr)
            offsets(n)=freq_epr(n)-freq_epr(centre);
        end
    end

end

% Consistency enforcement
function grumble(parameters,spin_system,B,euler_angles)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
    end
    if (~isnumeric(B))||(~isreal(B))||(~isscalar(B))
        error('B must be a real numeric scalar.');
    end
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||...
            (numel(euler_angles)~=3)
        error('euler_angles must be a three-element real vector.');
    end
end

