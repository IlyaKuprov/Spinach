%KEHL_NUCLEAR_CS_MATRIX Nuclear chemical-shift matrix for a Spinach spin.
%
%   Spinach architecture migration May 2026 Talos

function M=kehl_nuclear_cs_matrix(spin_system,spin_idx)
    M=(spin_system.inter.zeeman.matrix{spin_idx}/...
       spin_system.inter.basefrqs(spin_idx)-eye(3,3))*1e12;
end
