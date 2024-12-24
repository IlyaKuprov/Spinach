% Kronecker product of two matrices in a tensor train format. Syntax:
%
%                                c=kron(a,b)
%
% WARNING: the result is not the same as the flat matrix Kronecker pro-
%          duct (it is a row and column permutation away from it), but
%          the resulting order of elements is consistent with the out-
%          put of the tensor train vectorization (ttclass/vec) operati-
%          on output.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/kron.m>

function c=kron(a,b)

% Shrink a and b before going any further
a=shrink(a); b=shrink(b);

% Read sizes and ranks of the operands
[a_ncores,~]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a);
[b_ncores,~]=size(b.cores); b_ranks=ranks(b); b_sizes=sizes(b);

% Check consistency
if a_ncores~=b_ncores
    error('tensor train structures must be consistent.');
end

% Preallocate the result
new_cores=cell(a_ncores,1);

% Kron the cores
for nc=1:a_ncores
    core_of_a=reshape(a.cores{nc,1},[a_ranks(nc,1)*a_sizes(nc,1),a_sizes(nc,2)*a_ranks(nc+1,1)]);
    core_of_b=reshape(b.cores{nc,1},[b_ranks(nc,1)*b_sizes(nc,1),b_sizes(nc,2)*b_ranks(nc+1,1)]);
    core_of_c=reshape(kron(core_of_b,core_of_a),[a_ranks(nc,1),a_sizes(nc,1),b_ranks(nc,1),b_sizes(nc,1),a_sizes(nc,2),a_ranks(nc+1,1),b_sizes(nc,2),b_ranks(nc+1,1)]);
    new_cores{nc,1}=reshape(permute(core_of_c,[1 3 4 2 7 5 6 8]),[a_ranks(nc,1)*b_ranks(nc,1),b_sizes(nc,1)*a_sizes(nc,1),b_sizes(nc,2)*a_sizes(nc,2),a_ranks(nc+1,1)*b_ranks(nc+1,1)]);
end

% Write the output data structure
c=ttclass; c.coeff=b.coeff*a.coeff;
c.cores=reshape(new_cores,[a_ncores,1]);
if a.coeff==0 || b.coeff==0
    % The result is exactly zero
    c.tolerance=0;
else
    % The relative tolerances sum up
    c.tolerance=a.coeff*b.tolerance+b.coeff*a.tolerance;
end

end

% Gentlemen, you're trying to negotiate something you will never be able
% to negotiate. If negotiated, it will not be ratified. And if ratified,
% it will not work.
%
% Russell Bretherton, UK negotiator in Brussels, 1955

