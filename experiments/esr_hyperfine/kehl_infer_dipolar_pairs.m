%KEHL_INFER_DIPOLAR_PAIRS ENDOR nuclear pairs with non-zero couplings.
%
%   Spinach architecture migration May 2026 Talos

function pairs=kehl_infer_dipolar_pairs(spin_system,endor_spins)
    pairs=zeros(numel(endor_spins)*(numel(endor_spins)-1)/2,2);
    pair_count=0;
    for n=1:numel(endor_spins)
        for k=n+1:numel(endor_spins)
            spin_a=endor_spins(n);
            spin_b=endor_spins(k);
            if norm(kehl_coupling_matrix(spin_system,spin_a,spin_b),'fro')>0
                pair_count=pair_count+1;
                pairs(pair_count,:)=[spin_a,spin_b];
            end
        end
    end
    pairs=pairs(1:pair_count,:);
end
