% Bit-revert permutation for the tensor train operator. Syntax:
%
%                       tt=revert(tt)
% 
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/revert.m>

function tt=revert(tt)

% Read sizes and ranks
[ncores,ntrains]=size(tt.cores);

% Swap bond indices
for n=1:ntrains
    for k=1:ncores
        tt.cores{k,n}=permute(tt.cores{k,n}, [4,2,3,1]);
    end
end

% Revert the train direction
tt.cores=tt.cores(ncores:-1:1,:);

end

% Asking for efficiency and adaptability in the same program is 
% like asking for a beautiful and modest wife... we'll probably
% have to settle for one or the other.
%
% Gerald M. Weinberg, "The psychology of computer programming"

% #NHEAD #NGRUM