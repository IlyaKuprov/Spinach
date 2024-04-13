% Converts a tensor train representation of a matrix into a matrix.
% Syntax: 
%                           answer=full(ttrain)
%
% Note: the result can be huge, careless use would crash the system.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/full.m>

function answer=full(ttrain)

% Preallocate the result
answer=zeros(size(ttrain));

% Get object dimensions
[ncores,ntrains]=size(ttrain.cores);

% Get tensor ranks
ttranks=ranks(ttrain);

% Get mode sizes
modesizes=sizes(ttrain);

% Loop over the buffer
for n=1:ntrains
    
    % Multiply up the tensor train
    next_term=reshape(ttrain.cores{ncores,n},[ttranks(ncores,n),modesizes(ncores,1)*modesizes(ncores,2)]);
    for k=(ncores-1):(-1):1
        next_term=reshape(ttrain.cores{k,n},[ttranks(k,n)*modesizes(k,1)*modesizes(k,2),ttranks(k+1,n)])*next_term;
        next_term=reshape(next_term,[ttranks(k,n),modesizes(k,1),modesizes(k,2),prod(modesizes(k+1:ncores,1)),prod(modesizes(k+1:ncores,2))]);
        next_term=reshape(permute(next_term,[1 4 2 5 3]),[ttranks(k,n),prod(modesizes(k:ncores,1))*prod(modesizes(k:ncores,2))]);
    end
    next_term=reshape(next_term,[prod(modesizes(1:ncores,1)),prod(modesizes(1:ncores,2))]);
    
    % Add the matrix to the result
    answer=answer+ttrain.coeff(n)*next_term;
    
end

end

% We must conduct research and then accept the results. If they don't
% stand up to experimentation, Buddha's own words must be rejected.
%
% Dalai Lama XIV

