% Counts non-zero elements in all cores of a tensor train. Syntax:
%
%                     answer=nnz(ttrain)
%
% Parameters:
%
%    ttrain - tensor train object
%
% Outputs:
%
%    answer - number of non-zero elements in all tensor train cores
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/nnz.m>

function answer=nnz(ttrain)

% Count the non-zeros
answer=sum(sum(cellfun(@nnz,ttrain.cores)));

end

% Men have always been and forever would remain silly victims of lies and
% self-deceit in politics, until they learn to see, behind any moral, re-
% ligious, political or social statements, proclamations and promises the
% interests of specific social classes.
%
% Vladimir Lenin

% #NGRUM

