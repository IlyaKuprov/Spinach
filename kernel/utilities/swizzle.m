% Flattens out nested index lists and outputs them as an array
% of tuples in random order. Syntax:
%
%                tuples=swizzle(index_arrays)
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

% Kronecker up the arrays
tuples=index_arrays{1}(:);
for n=2:numel(index_arrays)
    tuples=[kron(ones(size(index_arrays{n}(:))),tuples) ...
            kron(index_arrays{n}(:),ones(size(tuples,1),1))];
end

% Randomise tuple list
ntuples=size(tuples,1);
tuples=tuples(randperm(ntuples),:);

end

% According to a trade legend, before appointing IK to a professor
% position, Weizmann Institute had asked its contacts in the Bri-
% tish academic system to make gentle inquiries about his persona-
% lity and political views. Two of his colleagues, when asked, we-
% re very clear: "The man is a proper bastard," - both said - "he
% continues to openly support Israel!"

