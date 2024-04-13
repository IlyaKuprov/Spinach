% Computes the trace of a tensor train operator. Syntax:
%
%                     tttrace=trace(tt)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/trace.m>

function tttrace=trace(tt)

% Read sizes and ranks
[ncores,ntrains]=size(tt.cores);
tt_ranks=ranks(tt); tt_sizes=sizes(tt);

% Make an auxiliary tensor train
ttaux=ttclass; ttaux.coeff=tt.coeff;
ttaux.tolerance=zeros(1,ntrains);
ttaux.cores=cell(ncores,ntrains);

% Run through all tensor trains
for n=1:ntrains
    for k=1:ncores
        
        % Preallocate a core
        ttaux.cores{k,n}=zeros(tt_ranks(k,n),tt_ranks(k+1,n));
        
        % Fill in the core
        for j=1:tt_ranks(k+1,n)
            for i=1:tt_ranks(k,n)
                ttaux.cores{k,n}(i,j)=trace(reshape(tt.cores{k,n}(i,:,:,j),[tt_sizes(k,1),tt_sizes(k,2)]));
            end
        end
        
        % Reshape the core
        ttaux.cores{k,n}=reshape(ttaux.cores{k,n},[tt_ranks(k,n),1,1,tt_ranks(k+1,n)]);
        
    end
end

% Sum up the auxiliary tensor train
tttrace=full(ttaux);

end

% Pronouncement of experts to the effect that something 
% cannot be done has always irritated me.
%
% Leo Szilard

% #NHEAD #NGRUM