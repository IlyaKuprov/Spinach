% Mean of elements of a tensor train representation of a
% matrix. Syntax:
%
%               answer=mean(ttrain,dim)
%
% The mean value is computed along the specified dimen-
% sion (dim=1 or dim=2).
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/mean.m>

function answer=mean(ttrain,dim)

% Get sizes and ranks
[ncores,ntrains]=size(ttrain.cores);
tt_ranks=ranks(ttrain);
tt_sizes=sizes(ttrain);

% If all dimensions are singleton, return a scalar immediately
if all(tt_sizes(:)==1), answer=full(ttrain); return; end

% In dim is omitted, choose first non-singleton dimension
% (this mimics the Matlab behaviour for matices) 
if nargin==1
    dim=0;
    for k=2:-1:1
        if ~all(tt_sizes(:,k)==1)
            dim=k;
        end
    end
end

% Make an auxiliary tensor train
answer=ttclass;
answer.coeff=ttrain.coeff;
answer.tolerance=zeros(1,ntrains);
answer.cores=cell(ncores,ntrains);

% Run through all tensor trains
for n=1:ntrains

    % Run through all cores
    for k=1:ncores
        
        % Sum the appropriate dimension in each core
        answer.cores{k,n}=sum(ttrain.cores{k,n},dim+1);
        
        % Reshape the core and divide the result by the correspondent dimension
        switch dim    
            case 1
                answer.cores{k,n}=reshape(answer.cores{k,n},[tt_ranks(k,n),1,tt_sizes(k,2),tt_ranks(k+1,n)])/tt_sizes(k,1);
            case 2
                answer.cores{k,n}=reshape(answer.cores{k,n},[tt_ranks(k,n),tt_sizes(k,1),1,tt_ranks(k+1,n)])/tt_sizes(k,2);
            otherwise
                error('incorrect dimension specificaton.');
        end
        
    end
    
end

% If all modes are singleton, transform to a scalar
if all(all(sizes(answer)==1)), answer=full(answer); end

end

% Fascism, nazism, communism, and socialism are only superficial 
% variations of the same monstrous theme - collectivism.
%
% Ayn Rand

