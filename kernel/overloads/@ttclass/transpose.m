% Transposes a tensor without complex conjugation. Syntax:
%
%                  ttrain=transpose(ttrain)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/transpose.m>

function ttrain=transpose(ttrain)

% Read tensor sizes and ranks 
[ncores,ntrains]=size(ttrain.cores);

% Swap the middle dimensions of all cores
for n=1:ntrains
    for k=1:ncores
        ttrain.cores{k,n}=permute(ttrain.cores{k,n},[1 3 2 4]);
    end
end

end

% "Public welfare" is the welfare of those who do not earn
% it; those who do, are entitled to no welfare.
%
% Ayn Rand, "Atlas Shrugged"

% #NHEAD #NGRUM