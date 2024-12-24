% Stretches arrays into vectors - useful for situations when the stand-
% ard (:) syntax is not available. Syntax:
%
%                              A=vec(A)
% 
% WARNING: for tensor trains this operation proceeds by stretching eve-
%          ry core of the tensor train. the result is not the same as
%          column-wise matrix stretching (it is an element permutation
%          away from it), but the resulting order of elements is consi-
%          stent with tensor train Kronecker product operation output.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/vec.m>

function A=vec(A)

% Decide how to proceed
if isa(A,'ttclass')
    
    % Read tensor train sizes and ranks
    [ncores,ntrains]=size(A.cores);
    ttm_ranks=ranks(A); ttm_sizes=sizes(A);
    
    % Reshape the cores
    for n=1:ntrains
        for k=1:ncores
            A.cores{k,n}=reshape(A.cores{k,n},[ttm_ranks(k,n),ttm_sizes(k,1)*ttm_sizes(k,2),1,ttm_ranks(k+1,n)]);
        end
    end
    
else
    
    % Use standard Matlab stretch
    A=reshape(A,[numel(A) 1]);
    
end
        
end

% Jean-Paul Sartre was sitting at a French cafe, revising his draft of Being
% and Nothingness. He said to the waitress, "I'd like a cup of coffee, plea-
% se, with no cream." The waitress replied, "I'm sorry, Monsieur, but we are
% out of cream. How about with no milk?"

