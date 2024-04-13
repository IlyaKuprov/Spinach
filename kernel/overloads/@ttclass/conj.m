% Conjugate the elements of the tensor train matrix. Syntax:
%
%                         tt=conj(tt)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/conj.m>

function tt=conj(tt)

% Read tensor train sizes and ranks
[ncores,ntrains]=size(tt.cores);

% Conjugate the cores
for n=1:ntrains
    for k=1:ncores
        tt.cores{k,n}=conj(tt.cores{k,n});
    end
end

% Conjugate the coefficients
tt.coeff=conj(tt.coeff);

end

% I asked God for a bike, but I know God 
% doesn't work that way. So I stole a bike
% and asked for forgiveness.
%
% Al Pacino

% #NGRUM #NHEAD