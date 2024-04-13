% Computes a Hermitian conjugate of a matrix in a tensor train
% representation. Syntax:
%
%                   ttrain=ctranspose(ttrain)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/ctranspose.m>

function ttrain=ctranspose(ttrain)

% Read tensor sizes and ranks 
[ncores,ntrains]=size(ttrain.cores);

% Swap the middle dimensions of all cores
for n=1:ntrains
    for k=1:ncores
        ttrain.cores{k,n}=permute(ttrain.cores{k,n},[1 3 2 4]);
    end
end

% Conjugate the result
ttrain=conj(ttrain);

end

% What gives the artist real prestige is his imitators.
%
% Igor Stravinsky

% #NGRUM #NHEAD