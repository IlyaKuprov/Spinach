% Wigner D function basis set and rotation generators required by
% the SLE module. Syntax:
%
%     [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)
%
% Parameters:
%
%     max_rank   - maximum L rank for Wigner D functions
%
%     int_ranks  - row vector of interaction ranks for which
%                  product superoperators are required; may
%                  be empty if only rotation generators are
%                  needed
%
% Outputs:
%
%     space_basis  -   lab space basis set descriptor, in
%                      [L M N] format, giving indices of
%                      each Wigner function in the basis.
%
%     Lx,Ly,Lz  - representations of lab space rotation
%                 generators in the Wigner function basis,
%                 to be used in the building of the lab
%                 space diffusion operator.
%
%     D         - a cell array with one element per interac-
%                 tion rank r in int_ranks, each a cell array
%                 of Wigner function product superoperators,
%                 corresponding to multiplication by D[r,M,N]
%                 of the basis Wigner functions, to be used
%                 in the building of the spin Hamiltonian
%                 operator; D{r} has dimensions (2r+1)x(2r+1)
%
% Automatic caching is implemented - the function would not re-
% compute operator sets that it can find on disk.
%
% Product superoperators at every interaction rank use the native
% MATLAB clebsch_gordan.m routine.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=sle_operators.m>

function [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)

% Check consistency
grumble(max_rank,int_ranks);

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-13));
cache_file=[own_path 'sle_operators_rank_' ...
            num2str(max_rank) '_int' ...
            sprintf('_%d',int_ranks) '.mat'];

% Check the cache
if exist(cache_file,'file')
    
    % Lift data from the cache if the file is already available
    load(cache_file,'space_basis','Lx','Ly','Lz','D');
    
else
    
    % Determine the spatial problem dimension
    basis_dim=(1+max_rank)*(1+2*max_rank)*(3+2*max_rank)/3;
    
    % Preallocate the spatial basis descriptor
    space_basis=zeros(basis_dim,3);
    
    % Populate the spatial basis descriptor
    for rank=0:max_rank
        
        % Determine the index extents
        extents_upper=rank*(2*rank-1)*(2*rank+1)/3+1;
        extents_lower=(1+rank)*(1+2*rank)*(3+2*rank)/3;
        
        % Assign Wigner function ranks
        space_basis(extents_upper:extents_lower,1)=rank;
        
        % Assign Wigner function left projection indices
        space_basis(extents_upper:extents_lower,2)=kron((rank:-1:-rank)',ones(2*rank+1,1));
        
        % Assign Wigner function right projection indices
        space_basis(extents_upper:extents_lower,3)=kron(ones(2*rank+1,1),(rank:-1:-rank)');
        
    end
    
    % Build spatial L+ generator: pick the states that can be raised
    source_states=space_basis(space_basis(:,2)<space_basis(:,1),:); new_dim=size(source_states,1);
    
    % Build spatial L+ generator: find out what they are raised into  
    destin_states=source_states+[zeros(new_dim,1) ones(new_dim,1) zeros(new_dim,1)];
    
    % Build spatial L+ generator: compute the corresponding matrix elements
    matrix_elements=sqrt(source_states(:,1).*(source_states(:,1)+1)-...
                         source_states(:,2).*(source_states(:,2)+1));
                     
    % Build spatial L+ generator: determine linear indices of source states
    L=source_states(:,1); M=source_states(:,2); N=source_states(:,3);
    source_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;
    
    % Build spatial L+ generator: determine linear indices of destination states
    L=destin_states(:,1); M=destin_states(:,2); N=destin_states(:,3);
    destin_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;
    
    % Build spatial L+ generator: form the L+ matrix
    Lp=sparse(destin_state_indices,source_state_indices,matrix_elements,basis_dim,basis_dim);
    
    % Build Lx and Ly rotation generators from L+
    Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;
    
    % Build Lz rotation generator
    Lz=spdiags(space_basis(:,2),0,basis_dim,basis_dim);
    
    % Do a clean-up
    clear('source_states','destin_states','matrix_elements','L','M','N',...
          'source_state_indices','destin_state_indices','Lp');

    % Preallocate product action matrices
    D=cell(1,max([int_ranks 0]));

    % Loop over the requested interaction ranks
    for int_rank=int_ranks

        % Preallocate the current rank
        n_prj=2*int_rank+1; D{int_rank}=cell(n_prj);

        % Build the source state list
        sources=kron(space_basis,ones(n_prj,1));

        % Loop over the left projection index
        for m=1:n_prj

            % Loop over the right projection index
            for n=1:n_prj

                % Get the indices of the acting Wigner function
                L1=int_rank; M1=int_rank+1-m; N1=int_rank+1-n;

                % Preallocate destination state list
                destinations=zeros(n_prj*basis_dim,3);

                % Loop over the basis set
                for k=1:basis_dim

                    % Get the source indices
                    L2=space_basis(k,1); M2=space_basis(k,2); N2=space_basis(k,3);

                    % Determine destination indices
                    ranks=((L2-L1):(L2+L1))';
                    prj_m=(M1+M2)*ones(size(ranks));
                    prj_n=(N1+N2)*ones(size(ranks));

                    % Build the destination state list
                    destinations((n_prj*(k-1)+1):(n_prj*k),:)=[ranks prj_m prj_n];

                end

                % Screen source and destination lists
                L=destinations(:,1); L2=sources(:,1);
                M=destinations(:,2); M2=sources(:,2);
                N=destinations(:,3); N2=sources(:,3);
                hit_list=(L>=max(0,abs(L2-L1)))&(abs(M)<=L)&(abs(N)<=L)&(L<=max_rank);
                L=L(hit_list); L1=L1*ones(size(L)); L2=L2(hit_list);
                M=M(hit_list); M1=M1*ones(size(M)); M2=M2(hit_list);
                N=N(hit_list); N1=N1*ones(size(N)); N2=N2(hit_list);

                % Scaling factor equivalent to normalizing the Wigner D functions
                normalization_factor=sqrt((2*L2+1)./(2*L+1));

                % Compute the structure coefficients of the Wigner D functions
                matrix_elements=zeros(size(L));
                parfor k=1:numel(L)
                    matrix_elements(k)=clebsch_gordan(L(k),M(k),L1(k),M1(k),L2(k),M2(k))*...
                                       clebsch_gordan(L(k),N(k),L1(k),N1(k),L2(k),N2(k))*...
                                       normalization_factor(k);
                end

                % Assign matrix elements
                destin_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;
                source_state_indices=L2.*(4*L2.^2+6*(L2-M2)+5)/3-M2-N2+1;
                D{int_rank}{m,n}=sparse(destin_state_indices,source_state_indices,...
                                        matrix_elements,basis_dim,basis_dim);

            end

        end

    end
        
    try % Try to save a cache record, but don't insist
        save(cache_file,'space_basis','Lx','Ly','Lz',...
                        'D','-v7.3'); drawnow;
    catch
        warning('Spinach installation appears to be write-protected');
    end
    
end
    
end

% Consistency enforcement
function grumble(max_rank,int_ranks)
if (numel(max_rank)~=1)||(~isnumeric(max_rank))||(~isreal(max_rank))||...
   (max_rank<1)||(mod(max_rank,1)~=0)
    error('max_rank must be a real positive integer.');
end
if (~isnumeric(int_ranks))||(~isreal(int_ranks))||...
   (~isempty(int_ranks)&&(~isrow(int_ranks)))||...
   any(mod(int_ranks,1)~=0)||any(int_ranks<1)||...
   (numel(unique(int_ranks))~=numel(int_ranks))
    error('int_ranks must be a row vector of distinct positive integers.');
end
end

% To invent, you need a good imagination and a pile of junk.
% 
% Thomas Edison


