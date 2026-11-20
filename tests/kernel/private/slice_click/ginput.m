% Supply one click in digitisation points and interrupt the next request.

function [x,y]=ginput(~)

% Interrupt the second request after the two traces have been plotted
persistent clicked
if ~isempty(clicked)
    error('test:clickComplete','One click completed.');
end

% Select different row and column indices to distinguish the traces
clicked=true;
x=3;
y=2;

end


