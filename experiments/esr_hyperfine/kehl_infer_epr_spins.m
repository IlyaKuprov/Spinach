%KEHL_INFER_EPR_SPINS Nuclei retained in the EPR orientation simulation.
%
%   Spinach architecture migration May 2026 Talos

function epr_spins=kehl_infer_epr_spins(spin_system,electron_idx,endor_spins)
    epr_spins=zeros(1,spin_system.comp.nspins);
    spin_count=0;
    for n=1:spin_system.comp.nspins
        if (n~=electron_idx) && (~ismember(n,endor_spins)) &&...
           (spin_system.comp.mults(n)>1)
            has_hfc=norm(kehl_coupling_matrix(spin_system,electron_idx,n),'fro')>0;
            has_nqi=norm(kehl_coupling_matrix(spin_system,n,n),'fro')>0;
            if has_hfc || has_nqi
                spin_count=spin_count+1;
                epr_spins(spin_count)=n;
            end
        end
    end
    epr_spins=epr_spins(1:spin_count);
end
