% Flattens out nested index lists and outputs them as an array
% of tuples in random order. This is useful for flattening nes-
% ted loops for parallel processing. Syntax:
%
%                 tuples=swizzle(index_arrays)
%
% Parameters:
%
%     index_arrays  - a cell array of row vectors 
%
% Outputs:
%
%     tuples        - a matrix of tuples in random or-
%                     der, with tuples listed as rows
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=swizzle.m>

function tuples=swizzle(index_arrays)

% Check consistency
grumble(index_arrays);

% Kronecker up the arrays
tuples=index_arrays{1}(:);
for n=2:numel(index_arrays)
    tuples=[kron(ones(size(index_arrays{n}(:))),tuples) ...
            kron(index_arrays{n}(:),ones(size(tuples,1),1))];
end

% Randomise the tuple list
ntuples=size(tuples,1); tuples=tuples(randperm(ntuples),:);

end

% Consistency enforcement
function grumble(index_arrays)
if ~iscell(index_arrays)
    error('index_arrays must be a cell array of row vectors.');
end
for n=1:numel(index_arrays)
    if (~isreal(index_arrays{n}))||(~isrow(index_arrays{n}))||...
       (any(mod(index_arrays{n},1)~=0,'all'))||...
       (any(index_arrays{n}<1,'all'))
        error('elements of index_arrays must be row vectors of positive integers.');
    end
end
end

% Acording to a conference rumour, before appointing IK to a tenured 
% full professor position, the Weizmann Institute had asked its peo-
% ple in the UK academic system to make inquiries about his persona-
% lity and political views. Two of IK's former colleagues at Oxford,
% without realising who's asking, had been very clear: "The man is a
% proper bastard," - both said - "for goodness sake, he continues to 
% openly support Israel!"

