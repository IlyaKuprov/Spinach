% Returns a unit object of the same type as whatever is supplied.
%
%                         A=unit_like(A)
%
% Full square matrices, sparse matrices and tensor train represen-
% tations of square matrices are supported.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/unit_like.m>

function A=unit_like(A)

if isa(A,'ttclass')
    
    % Unit tensor train of the same topology
    mode_sizes=sizes(A);
    if all(mode_sizes(:,1)==mode_sizes(:,2))
        core=cell(A.ncores,1);
        for k=1:A.ncores
            core{k}=eye(mode_sizes(k,1));
        end
        A=ttclass(1,core,0);
    else
        error('tensor train does not represent a square matrix.');
    end
    
elseif ismatrix(A)&&(size(A,1)==size(A,2))&&issparse(A)
    
    % Unit sparse matrix of the same dimension
    A=speye(size(A));
    
elseif ismatrix(A)&&(size(A,1)==size(A,2))&&(~issparse(A))
    
    % Unit dense matrix of the same dimension
    A=eye(size(A));
    
else
    
    % Complain and bomb out
    error('the input is not a square matrix or a representation thereof.');
    
end

end

% Briefly stated, the Gell-Mann Amnesia effect is as follows. You open the
% newspaper to an article on some subject you know well. You read the arti-
% cle and see the journalist has absolutely no understanding of either the
% facts or the issues. Often, the article is so wrong it actually presents
% the story backward -- reversing cause and effect. In any case, you read
% with exasperation or amusement the multiple errors in a story, and then
% turn the page to national or international affairs, and read as if the
% rest of the newspaper was somehow more accurate about Palestine than the
% baloney you just read. You turn the page, and forget what you know.
%
% Michael Crichton

