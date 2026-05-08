% ENDOR nuclear pairs with non-zero couplings. Syntax:
%
%      pairs=kehl_infer_dipolar_pairs(spin_system,endor_spins)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   endor_spins      - Spinach indices of ENDOR nuclei.
%
% Outputs:
%
%   pairs            - two-column array of coupled ENDOR spin pairs.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_infer_dipolar_pairs.m>

function pairs=kehl_infer_dipolar_pairs(spin_system,endor_spins)

    % Check consistency
    grumble(spin_system,endor_spins);

    % Preallocate the largest possible pair list
    pairs=zeros(numel(endor_spins)*(numel(endor_spins)-1)/2,2);
    pair_count=0;

    % Keep nuclear pairs with non-zero Spinach couplings
    for n=1:numel(endor_spins)
        for k=(n+1):numel(endor_spins)
            spin_a=endor_spins(n);
            spin_b=endor_spins(k);
            if norm(kehl_coupling_matrix(spin_system,spin_a,spin_b),'fro')>0
                pair_count=pair_count+1;
                pairs(pair_count,:)=[spin_a,spin_b];
            end
        end
    end

    % Trim unused rows
    pairs=pairs(1:pair_count,:);

end

% Consistency enforcement
function grumble(spin_system,endor_spins)
    if (~isstruct(spin_system))||(~isfield(spin_system,'inter'))||...
            (~isfield(spin_system.inter,'coupling'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(endor_spins))||(~isvector(endor_spins))||...
            any(endor_spins<1)||any(mod(endor_spins,1)~=0)
        error('endor_spins must be a vector of positive integers.');
    end
end

