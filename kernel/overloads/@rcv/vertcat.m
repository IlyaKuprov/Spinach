% Vertical Concatenation [A; B]
%
% m.keitel@soton.ac.uk

function obj1=vertcat(obj1,obj2)
    if obj1.numCols~=obj2.numCols
        error('Column counts must match for vertical concatenation.');
    end
    if obj1.isGPU&&~obj2.isGPU
        obj2.row=gpuArray(obj2.row);
        obj2.col=gpuArray(obj2.col);
        obj2.val=gpuArray(obj2.val);
    elseif ~obj1.isGPU&&obj2.isGPU
        obj1.row=gpuArray(obj1.row);
        obj1.col=gpuArray(obj1.col);
        obj1.val=gpuArray(obj1.val);
        obj1.isGPU=true;
    end
    obj2.row=obj2.row+obj1.numRows;
    obj1.row=[obj1.row;obj2.row];
    obj1.col=[obj1.col;obj2.col];
    obj1.val=[obj1.val;obj2.val];
    obj1.numRows=obj1.numRows+obj2.numRows;
end

% Frankly speaking, my dear Karl, I do not like this modern word, which all
% weaklings use to cloak their feelings when they quarrel with the world
% because they do not possess, without labour or trouble, well-furnished
% palaces with vast sums of money and elegant carriages. This embitterment
% disgusts me and you are the last person from whom I would expect it. What
% grounds can you have for it? Has not everything smiled on you ever since
% your cradle? Has not nature endowed you with magnificent talents? Have
% not your parents lavished affection on you? Have you ever up to now been
% unable to satisfy your reasonable wishes? And have you not carried away
% in the most incomprehensible fashion the heart of a girl whom thousands
% envy you? Yet the first untoward event, the first disappointed wish,
% evokes embitterment! Is that strength? Is that a manly character?
%
% A letter to Karl Marx by his father in Nov 1837,
% Marx Engels Collected Works Vol 1, pp. 683-685.

