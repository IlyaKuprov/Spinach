% Mimics the diag behavior for tensor train matrix. Syntax:
%
%                            tt=diag(tt)
%
% If tt is a square matrix, return a vector by computing diag of every core.
% If tt is a vector (one mode size is ones), return a diagonal matrix.
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/diag.m>

function tt=diag(tt)

% Read tensor train sizes and ranks
[ncores,ntrains]=size(tt.cores);
r=ranks(tt); sz=sizes(tt);

% Decide the dimensions
if all(sz(:,1)==ones(ncores,1)) || all(sz(:,2)==ones(ncores,1))
    
    % Vector on input, diagonal matrix on output
    for n=1:ntrains
        for k=1:ncores
            vector_core=tt.cores{k,n};
            matrix_size=max(sz(k,1),sz(k,2));
            matrix_core=zeros(r(k,n),matrix_size,matrix_size,r(k+1,n));
            for p=1:r(k,n)
                for q=1:r(k+1,n)
                    matrix_core(p,:,:,q)=diag(reshape(vector_core(p,:,:,q),[sz(k,1),sz(k,2)]));
                end
            end
            tt.cores{k,n}=matrix_core;
        end
    end
    
else
    
    % Matrix on input, column vector on output
    if all(sz(:,1)==sz(:,2))
        for n=1:ntrains
            for k=1:ncores
                matrix_core=tt.cores{k,n};
                vector_size=min(sz(k,1),sz(k,2));
                vector_core=zeros(r(k,n),vector_size,1,r(k+1,n));
                for p=1:r(k,n)
                    for q=1:r(k+1,n)
                        vector_core(p,:,1,q)=diag(reshape(matrix_core(p,:,:,q),[sz(k,1),sz(k,2)]));
                    end
                end
                tt.cores{k,n}=vector_core;
            end
        end
    else
        error('Input should be either a square matrix or a vector.');
    end
    
end

end

% The only mistake [the famous criminal finacier] Bernie Madoff 
% made was to promise returns in this life.
%
% Harry Kroto

