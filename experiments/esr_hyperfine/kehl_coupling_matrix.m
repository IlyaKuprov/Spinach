%KEHL_COUPLING_MATRIX Spinach coupling matrix converted to Hz.
%
%   Spinach architecture migration May 2026 Talos

function M=kehl_kehl_coupling_matrix(spin_system,row,col)
    M=kehl_matrix_from_cell(spin_system.inter.coupling.matrix,row,col)/(2*pi);
end
