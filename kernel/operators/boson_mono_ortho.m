function B = boson_mono_ortho(nlevels)

% Calling the bosonic monomials to be orthogonalized
Bm = boson_mono(nlevels);

% Extracting number of elements
N = numel(Bm); 

% Preallocating cell array for orthogonalized bosonic monomials
B = cell(N,1);

% Tolerance for numerical linear dependence
tol = 1e-12;

% Gram-Schmidt orthogonalization
for n = 1:N
    total_proj = sparse(zeros(nlevels)); 
    for k = 1:(n-1)
        denom = hdot(B{k},B{k});
        % Checking if the norm in the denominator is stable
        if abs(denom)<tol
            continue;
        end    
        sub_proj = (hdot(B{k},Bm{n})/denom)*B{k};
        total_proj = total_proj + sub_proj;
    end
        B{n} = Bm{n} - total_proj; % Storing orthogonal monomials 

end

end