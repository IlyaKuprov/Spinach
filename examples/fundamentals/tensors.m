% Test the conversion from Stevens operator coefficients to
% irreducible spherical tensor coefficients.
%
% i.kuprov@soton.ac.uk

function tensors()

% Test up to rank 6 on spin 15/2
max_rank=6; mult=16;

% Generate a random set of Stevens operator coefficients
for k=1:max_rank, r{k}=rand(2*k+1,1)-0.5; end %#ok<AGROW>

% Build the linear combination
lin_comb_a=sparse(0);
for k=1:max_rank
    for q=1:(2*k+1)
        lin_comb_a=lin_comb_a+r{k}(q)*stevens(mult,k,q-k-1);
    end
end
    
% Translate the coefficients into ISTs
for k=1:max_rank, r{k}=stev2sph(k,r{k}); end
   
% Build the linear combination
lin_comb_b=sparse(0);
for k=1:max_rank
    T=irr_sph_ten(mult,k);
    for q=1:(2*k+1)
        lin_comb_b=lin_comb_b+r{k}(q)*T{q};
    end
end

% Subtract the matrices and check the norm
if norm(lin_comb_b-lin_comb_a,1)<1e-6
    disp('Test passed.');
else
    disp('Test FAILED.');
    disp(full(lin_comb_a));
    disp(full(lin_comb_b));
end

end

