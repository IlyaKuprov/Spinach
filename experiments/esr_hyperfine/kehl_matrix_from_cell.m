% Matrix lookup from Spinach interaction cell arrays. Syntax:
%
%      M=kehl_matrix_from_cell(cells,row,col)
%
% Parameters:
%
%   cells            - Spinach interaction cell array.
%   row,col          - spin indices to be queried.
%
% Outputs:
%
%   M                - interaction matrix, or a zero matrix when absent.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_matrix_from_cell.m>

function M=kehl_matrix_from_cell(cells,row,col)

    % Check consistency
    grumble(cells,row,col);

    % Return a zero matrix for absent interactions
    M=zeros(3,3);
    if isempty(cells)
        return
    end

    % Prefer direct storage, then try the transposed cell
    if size(cells,1)>=row&&size(cells,2)>=col&&~isempty(cells{row,col})
        M=cells{row,col};
    elseif size(cells,1)>=col&&size(cells,2)>=row&&~isempty(cells{col,row})
        M=cells{col,row};
    end

end

% Consistency enforcement
function grumble(cells,row,col)
    if ~iscell(cells)
        error('cells must be a cell array.');
    end
    if (~isnumeric(row))||(~isscalar(row))||(row<1)||(mod(row,1)~=0)
        error('row must be a positive integer.');
    end
    if (~isnumeric(col))||(~isscalar(col))||(col<1)||(mod(col,1)~=0)
        error('col must be a positive integer.');
    end
end

