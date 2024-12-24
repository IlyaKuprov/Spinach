% Returns the bond dimensions of a tensor train. Syntax:
%
%                    ttranks=ranks(ttrain)
%
% The output is (ncores+1) by (ntrains) array. The first
% and the last elements for each train are 1.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/ranks.m>

function ttranks=ranks(ttrain)

% Get core array dimensions
[ncores,ntrains]=size(ttrain.cores);

% Preallocate the answer
ttranks=zeros(ncores+1,ntrains);

% Loop over the buffer
for n=1:ntrains
    
    % Extract the ranks
    for k=1:ncores
        ttranks(k,n)=size(ttrain.cores{k,n},1);
    end
    ttranks(ncores+1,n)=size(ttrain.cores{ncores,n},4);
    
end

end

% I refrain from publishing for fear that disputes and controversies
% may be raised against me by ignoramuses.
%
% Isaac Newton, in a letter to Leibniz

