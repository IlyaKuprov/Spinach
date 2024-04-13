% Haramard dot product between two tensor train matrices. Syntax:
%
%                         c=hdot(a,b)
%
% Parameters:
%
%    a,b  - tensor train objects representing numerical
%           arrays of the same dimensions and having
%           the same internal topology
%
% Outputs:
%
%    c    - Hadamard product of a and b, a scalar
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/hdot.m>

function c=hdot(a,b)

% Check consistency
grumble(a,b);

% Read topology and initialize the answer
[ncores,ntrains_a]=size(a.cores); ranks_a=ranks(a);
[~     ,ntrains_b]=size(b.cores); ranks_b=ranks(b);
mode_sizes=sizes(a); c=0;

% Loop over TT buffers
for na=1:ntrains_a
    for nb=1:ntrains_b
        
        % Multiply coefficients
        x=conj(a.coeff(na))*b.coeff(nb);
        
        % Loop over TT cores and compute dot product 
        for k=1:ncores
            core_b=reshape(b.cores{k,nb},[ranks_b(k,nb),mode_sizes(k,1)*mode_sizes(k,2)*ranks_b(k+1,nb)]);
            core_b=x*core_b;
            core_b=reshape(core_b,[ranks_a(k,nb)*mode_sizes(k,1)*mode_sizes(k,2),ranks_b(k+1,nb)]);      
            core_a=reshape(a.cores{k,na},[ranks_a(k,na)*mode_sizes(k,1)*mode_sizes(k,2),ranks_a(k+1,na)]);
            x=core_a'*core_b;          
        end
        
        % Add to the total
        c=c+x;
        
    end
end

end

% Consistency enforcement
function grumble(a,b)
if ~isa(a,'ttclass') || ~isa(b,'ttclass') 
    error('both arguments should be tensor trains.')
end
if size(a.cores,1)~=size(b.cores,1)
    error('the number of cores in the arguments must be the same.')
end
if ~all(all(sizes(a)==sizes(b)))
    error('tensor train mode sizes are not consistent.')
end
end

% The higher we soar, the smaller we appear to 
% those who cannot fly.
%
% Friedrich Nietzsche,
% "Thus Spoke Zarathustra"

