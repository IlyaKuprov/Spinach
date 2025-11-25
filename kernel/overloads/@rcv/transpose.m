% The transpose of an RCV sparse matrix. Syntax:
%
%                 A=transpose(A)
%
% Parameters:
%
%    A   - an RCV sparse matrix
%
% Outputs:
%
%    A   - transposed RCV matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/transpose.m>

function A=transpose(A)

% Check consistency
grumble(A);

% Efficiently swap rows and columns
[A.col,A.row]=deal(A.row,A.col);

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
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

