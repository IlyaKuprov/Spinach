% Magnetic-quantum-number basis table for EPR nuclei. Syntax:
%
%      basis_table=kehl_ori_field_build_space(spin_q)
%
% Parameters:
%
%   spin_q           - spin quantum numbers.
%
% Outputs:
%
%   basis_table      - basis vectors in columns.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ori_field_build_space.m>

function basis_table=kehl_ori_field_build_space(spin_q)

    % Check consistency
    grumble(spin_q);

    % Initialise the empty basis table
    basis_table=[];

    % Initialise the dimension counter
    dim=1;

    % Append one magnetic quantum number column at a time
    for n=1:length(spin_q)
        basis_block=[];
        for k=-spin_q(n):spin_q(n)
            basis_block=[basis_block; basis_table k*ones(dim,1)];
        end
        basis_table=basis_block;
        dim=size(basis_table,1);
    end

end

% Consistency enforcement
function grumble(spin_q)
    if ~isnumeric(spin_q)
        error('spin_q must be numeric.');
    end
end

