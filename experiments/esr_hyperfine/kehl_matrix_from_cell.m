%KEHL_MATRIX_FROM_CELL Matrix lookup from Spinach interaction cell arrays.
%
%   Spinach architecture migration May 2026 Talos

function M=kehl_kehl_matrix_from_cell(cells,row,col)
    M=zeros(3,3);
    if isempty(cells)
        return
    end
    if size(cells,1)>=row && size(cells,2)>=col && ~isempty(cells{row,col})
        M=cells{row,col};
    elseif size(cells,1)>=col && size(cells,2)>=row && ~isempty(cells{col,row})
        M=cells{col,row};
    end
end
