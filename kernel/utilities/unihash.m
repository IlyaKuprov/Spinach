% Hash table based stable duplicate row eliminator,
% for use with large sparse matrices where Matlab's
% unique(...,'rows') is too slow. Syntax:
%
%                   A=unihash(A)
%
% Parameters:
%
%      A - a large and sparse matrix
%
% Outputs:
%
%      A - same matrix with duplicate 
%          rows deleted, keeping the
%          first occurrence of each
%
% ilya.kuprov@weizmann.ac.il
% 
% <https://spindynamics.org/wiki/index.php?title=unihash.m>

function A=unihash(A)

% Build an MD5 hash table
hash_table=repmat(' ',[size(A,1) 32]);
parfor k=1:size(A,1)
    hash_table(k,:)=md5_hash(A(k,:));
end

% Redundant row index using a hash table
[~,idx]=unique(hash_table,'rows','stable');

% Elimination
A=A(idx,:);

end

% Вот полно статей: норм введение, норм постановка задачи. А потом 
% начинается ХУЙНЯ. Лютая, бестолковая, бесполезная ХУЙНЯ. И такие
% мол вот: мы тут внесли, расширили, исследовали. Солнышки -- нет!
% Вы ХУЙНЮ нагородили по данному вопросу. Лютую, бестолковую, бес-
% полезную ХУЙНЮ. Теоретики особенно этим отличаются.
%
% Anonymous post in
% an academic forum

