%KEHL_ELECTRON_G_MATRIX Electron g matrix for a Spinach spin.
%
%   Spinach architecture migration May 2026 Talos

function M=kehl_electron_g_matrix(spin_system,spin_idx)
    M=spin_system.inter.zeeman.matrix{spin_idx}*...
      spin_system.tols.freeg/spin_system.inter.basefrqs(spin_idx);
end
