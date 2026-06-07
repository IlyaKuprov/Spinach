% Returns mode sizes (physical dimensions of each core) of
% a tensor train. Syntax:
%
%                   modesizes=sizes(tt)
%
% Parameters:
%
%    tt - tensor train object
%
% Outputs:
%
%    modesizes - ncores by 2 array of physical dimensions
%                of tensor train cores
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/sizes.m>

function modesizes=sizes(tt)

% Determine the number of cores
ncores=size(tt.cores,1);

% Preallocate the answer
modesizes=zeros(ncores,2);

% Fill in the answer
for k=1:ncores
    modesizes(k,1)=size(tt.cores{k,1},2);
    modesizes(k,2)=size(tt.cores{k,1},3);
end

end

% Computer models are no different from fashion models: seductive,
% unreliable, easily corrupted, and they lead sensible people to
% make fools of themselves.
% 
% Jim Hacker, in "Yes, Prime Minister" (a BBC documentary)

% #NGRUM

