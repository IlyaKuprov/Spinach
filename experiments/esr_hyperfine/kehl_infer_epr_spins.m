% Nuclei retained in the EPR orientation simulation. Syntax:
%
%      epr_spins=kehl_infer_epr_spins(spin_system,electron_idx,endor_spins)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   electron_idx     - electron spin index.
%   endor_spins      - Spinach indices of ENDOR nuclei.
%
% Outputs:
%
%   epr_spins        - indices of auxiliary EPR nuclei.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_infer_epr_spins.m>

function epr_spins=kehl_infer_epr_spins(spin_system,electron_idx,endor_spins)

    % Check consistency
    grumble(spin_system,electron_idx,endor_spins);

    % Preallocate the largest possible EPR-nucleus list
    epr_spins=zeros(1,spin_system.comp.nspins);
    spin_count=0;

    % Keep nuclei with hyperfine, or quadrupolar couplings
    for n=1:spin_system.comp.nspins
        if (n~=electron_idx)&&(~ismember(n,endor_spins))&&...
                (spin_system.comp.mults(n)>1)
            has_hfc=norm(kehl_coupling_matrix(spin_system,electron_idx,n),'fro')>0;
            has_nqi=norm(kehl_coupling_matrix(spin_system,n,n),'fro')>0;
            if has_hfc||has_nqi
                spin_count=spin_count+1;
                epr_spins(spin_count)=n;
            end
        end
    end

    % Trim unused entries
    epr_spins=epr_spins(1:spin_count);

end

% Consistency enforcement
function grumble(spin_system,electron_idx,endor_spins)
    if (~isstruct(spin_system))||(~isfield(spin_system,'comp'))||...
            (~isfield(spin_system.comp,'nspins'))||(~isfield(spin_system.comp,'mults'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(electron_idx))||(~isscalar(electron_idx))||...
            (electron_idx<1)||(mod(electron_idx,1)~=0)
        error('electron_idx must be a positive integer.');
    end
    if (~isnumeric(endor_spins))||(~isvector(endor_spins))||...
            any(endor_spins<1)||any(mod(endor_spins,1)~=0)
        error('endor_spins must be a vector of positive integers.');
    end
end

