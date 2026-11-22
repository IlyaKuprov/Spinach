% Converts occupied cells of a pair table to ordered records. Syntax:
%
%                     pairs=cell_records(table)
%
% Parameters:
%
%    table - two-dimensional cell array
%
% Outputs:
%
%    pairs - cell rows {from,to,value} in column-major order
%
% talos@spindynamics.org

function pairs=cell_records(table)

% Check consistency
grumble(table);

% Visit one column at a time to avoid a square occupancy mask
columns=cell(size(table,2),1);
for k=1:size(table,2)
    column=table(:,k);
    rows=find(~cellfun(@isempty,column));
    columns{k}=[num2cell(rows(:)) repmat({k},numel(rows),1) reshape(column(rows),[],1)];
end
pairs=vertcat(cell(0,3),columns{:});

end

% Consistency enforcement
function grumble(table)
if (~iscell(table))||(~ismatrix(table))
    error('table must be a two-dimensional cell array.');
end
end


