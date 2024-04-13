% Performs right-to-left SVD recompression for a tensor train. This
% should not be called directly, use shrink.m instead. Syntax:
%
%                        ttout=truncate(tt)
%
% Parameters:
%
%     tt - a tensor train object with tt.ntrains=1 
%          and orthogonalised left-to-right
%
% Outputs:
%
%     ttout - a tensor train object, orthogonalised 
%             right-to-left
%
% Note: approximation tolerance (in Frobenius norm) is read from
%       tt.tolerance property.
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/truncate.m>

function ttout=truncate(tt)

% Check consistency
grumble(tt);

% Read tensor ranks and dimensions
sz=tt.sizes; rnk=tt.ranks; d=tt.ncores;

% Preallocate the result
ttout=tt; r=rnk(:,1);

% Define the relative approximation accuracy
eps=tt.tolerance(1,1)/tt.coeff(1,1);
eps=eps/sqrt(d);

% SVD cores right-to-left
for k=d:-1:2
    C=reshape(ttout.cores{k,1}, [r(k), sz(k,1)*sz(k,2)*r(k+1)]);
    B=reshape(ttout.cores{k-1,1}, [r(k-1)*sz(k-1,1)*sz(k-1,2), r(k)]);
    [U,S,V]=svd(C,'econ'); S=diag(S);
    rnew=frob_chop(S,eps*norm(S,'fro'));
    U=U(:,1:rnew); S=S(1:rnew); V=V(:,1:rnew);
    ttout.cores{k,1}=reshape(V', [rnew, sz(k,1), sz(k,2), r(k+1)]);
    ttout.cores{k-1,1}=reshape(B*U*diag(S), [r(k-1), sz(k-1,1), sz(k-1,2), rnew]);
    r(k)=rnew;
end

% Normalise the output
nrm=norm(ttout.cores{1,1}(:),2);
ttout.cores{1,1}=ttout.cores{1,1}/nrm;
ttout.coeff(1,1)=ttout.coeff(1,1)*nrm;

end

% Consistency enforcement
function grumble(tt)
if tt.ntrains>1
    error('Please sum all tensor trains before calling truncate');
end
end

% "Even with basic 128-bit hash functions, such as MD5, [a hash 
% collision] is a vanishingly rare event in the physical scien-
% ces context: a useful rule of thumb in modern physics is that
% any probability smaller than that of the researcher committ-
% ing suicide (approximately 9.7e-5 per year in the UK in 2013)
% is negligible."
%
% http://dx.doi.org/10.1063/1.4949534

