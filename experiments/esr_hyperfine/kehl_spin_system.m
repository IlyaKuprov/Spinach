% Reduced Spinach spin system for Kehl ENDOR kernels. Syntax:
%
%      spin_system=kehl_spin_system(isotopes,n_spin_systems)
%
% Parameters:
%
%   isotopes         - electron isotope followed by ENDOR nuclear isotopes.
%   n_spin_systems   - number of separated spin systems.
%
% Outputs:
%
%   spin_system      - Spinach spin system object in Zeeman-Hilbert formalism.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_spin_system.m>

function spin_system=kehl_spin_system(isotopes,n_spin_systems)

    % Default is a single compact ENDOR spin system
    if nargin<2
        n_spin_systems=1;
    end

    % Check consistency
    grumble(isotopes,n_spin_systems);

    % Preserve the separated-nucleus calculation dimension
    if n_spin_systems>1
        isotopes=isotopes(1:2);
    end

    % Set quiet local Spinach constructor input
    sys.magnet=1;
    sys.isotopes=isotopes;
    sys.output='hush';
    sys.disable={'hygiene'};

    % No interactions are needed for field-independent operator generation
    inter=[];

    % Use the Hilbert-space Zeeman basis expected by the Kehl kernels
    bas.formalism='zeeman-hilb';
    bas.approximation='none';

    % Build the reduced Spinach spin system
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

end

% Consistency enforcement
function grumble(isotopes,n_spin_systems)
    if (~iscell(isotopes))||isempty(isotopes)
        error('isotopes must be a non-empty cell array.');
    end
    if any(~cellfun(@ischar,isotopes))
        error('isotopes must contain character strings.');
    end
    if (~isnumeric(n_spin_systems))||(~isscalar(n_spin_systems))||...
            (n_spin_systems<1)||mod(n_spin_systems,1)~=0
        error('n_spin_systems must be a positive integer.');
    end
    if (n_spin_systems>1)&&(numel(isotopes)<2)
        error('at least one ENDOR nucleus is required for separated spin systems.');
    end
end

