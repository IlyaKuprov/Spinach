% Sum of elements of a tensor train representation of a matrix. Syntax:
%
%                        answer=sum(ttrain,dim)
%
% The sum is taken along the specified dimension (dim=1 or dim=2).
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/sum.m>

function answer=sum(ttrain,dim)

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
        
        % Reshape the core
        switch dim    
            case 1
                answer.cores{k,n}=reshape(answer.cores{k,n},[tt_ranks(k,n),1,tt_sizes(k,2),tt_ranks(k+1,n)]);
            case 2
                answer.cores{k,n}=reshape(answer.cores{k,n},[tt_ranks(k,n),tt_sizes(k,1),1,tt_ranks(k+1,n)]);
            otherwise
                error('incorrect dimension specificaton.');
        end
        
    end
    
end

% If all modes are singleton, transform to a scalar
if all(all(sizes(answer)==1)), answer=full(answer); end

end

% The sledge rests in summer, the cart rests in 
% winter, but the horse never rests.
%
% A Russian saying

