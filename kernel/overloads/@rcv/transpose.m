% The transpose of an RCV sparse matrix. Syntax:
%
%                 obj=transpose(obj)
%
% Parameters:
%
%    obj   - an RCV sparse matrix
%
% Outputs:
%
%    obj   - transposed RCV matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/transpose.m>

function obj=transpose(obj)

% Check consistency
grumble(obj);

% Efficiently swap rows and columns
[obj.col,obj.row]=deal(obj.row,obj.col);

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Я Шойгу. Значит, объясняю. Если вы такой хороший хозяин, что у вас
% котёнок умудрился свалиться в мусоропровод, то, во-первых, не надо
% прыгать и вопить "Барсик, милый, сука, держись!" Потому что держаться там
% не за что. Не надо пытаться пробить мусоропровод кувалдой, глухой
% котёнок-идиот не будет за это вам благодарен. Не нужно пытаться достать
% его проволочным крюком, дырявый котёнок-инвалид тоже не будет вам лучшим
% другом. Нужно просто взять в дрожащую руку телефон, а недрожащей набрать
% 911. Я приеду, разверну котёнкодоставатель и все сделаю за минуту.
%
% Attributed to Sergei Shoigu during his time as
% the Head of the Russian Emergencies Command

