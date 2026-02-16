function B_new = boson_mono_ortho(nlevels)

% Calling the bosonic monomials to be orthogonalized
B_old = boson_mono(nlevels);

% Extracting number of elements
N = numel(B_old); 

% Preallocating cell array for orthogonalized bosonic monomials
B_new = cell(N,1);

% Gram-Schmidt orthogonalization
for n = 1:N
    total_proj = spalloc(nlevels,nlevels,0); 
    for k = 1:(n-1)
        denom = hdot(B_new{k},B_new{k});
        sub_proj = (hdot(B_new{k},B_old{n})/denom)*B_new{k};
        total_proj = total_proj + sub_proj;
    end
        B_new{n} = B_old{n} - total_proj; % Storing orthogonal monomials 

end

end