% Wigner D function basis set and rotation generators required by 
% the SLE module. Syntax:
%
%        [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank)
%
% Input parameters:
%
%     max_rank   - maximum L rank for Wigner D functions
%
% Output parameters:
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
%     D         - a cell array of Wigner function product
%                 superoperators, corresponding to multi-
%                 plication by D[2,M,N] of the basis Wig-
%                 ner functions, to be used in the build-
%                 ing of the spin Hamiltonian operator.
%
% Automatic caching is implemented - the function would not re-
% compute operator sets that it can find on disk.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=sle_operators.m>

function [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank)

% Check consistency
grumble(max_rank);

% Generate cache record name
own_path=mfilename('fullpath'); 
own_path=own_path(1:(end-13));
cache_file=[own_path 'sle_operators_rank_' ...
            num2str(max_rank) '.mat'];

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

    % Preallocate product action matrices for second-rank Wigner functions
    D=cell(5);

    % Build the source state list
    sources=kron(space_basis,ones(5,1));
    
    % Loop over the left projection index
    for m=1:5
        
        % Loop over the right projection index
        for n=1:5
            
            % Get the indices of the acting Wigner function
            L1=2; M1=3-m; N1=3-n;
            
            % Preallocate destination state list
            destinations=zeros(5*basis_dim,3);
            
            % Loop over the basis set
            for k=1:basis_dim
                
                % Get the source indices
                L2=space_basis(k,1); M2=space_basis(k,2); N2=space_basis(k,3);
                
                % Determine destination indices
                ranks=((L2-L1):(L2+L1))';
                prj_m=(M1+M2)*ones(size(ranks));
                prj_n=(N1+N2)*ones(size(ranks));
                
                % Build the destination state list
                destinations((5*(k-1)+1):(5*(k-1)+5),:)=[ranks prj_m prj_n];
                
            end
            
            % Screen source and destination lists
            L=destinations(:,1); L2=sources(:,1); 
            M=destinations(:,2); M2=sources(:,2);
            N=destinations(:,3); N2=sources(:,3);
            hit_list=(L>=max(0,abs(L2-L1)))&(abs(M)<=L)&(abs(N)<=L)&(L<=max_rank);
            L=L(hit_list); L1=L1*ones(size(L)); L2=L2(hit_list);
            M=M(hit_list); M1=M1*ones(size(M)); M2=M2(hit_list);
            N=N(hit_list); N1=N1*ones(size(N)); N2=N2(hit_list);
            
            % Scaling factor for the structure coefficients that is equivalent 
            % to normalizing the Wigner D functions used
            normalization_factor=sqrt((2*L2+1)./(2*L+1));
            
            % Compute the structure coefficients of the Wigner D functions
            matrix_elements=clebsch_gordan_bypass(L,M,L1,M1,L2,M2).*...
                            clebsch_gordan_bypass(L,N,L1,N1,L2,N2).*...
                            normalization_factor;
                
            % Assign matrix elements
            destin_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;
            source_state_indices=L2.*(4*L2.^2+6*(L2-M2)+5)/3-M2-N2+1;
            D{m,n}=sparse(destin_state_indices,source_state_indices,matrix_elements,basis_dim,basis_dim);
                
        end
        
    end
        
    try % Try to save a cache record, but don't insist
        save(cache_file,'space_basis','Lx','Ly','Lz','D','-v7.3');
    catch
        warning('Spinach installation appears to be write-protected');
    end
    
end
    
end

% Hard-coded Clebsch-Gordan formulae for the L1=2 case
function cg=clebsch_gordan_bypass(L_array,~,~,M1_array,L2_array,M2_array)

% Preallocate the answer
cg=zeros(size(L_array));

parfor n=1:numel(L_array)
    
    % Get local variables
    L=L_array(n); M1=M1_array(n); L2=L2_array(n); M2=M2_array(n);

    % Enumerate all cases explicitly
    switch M1
        case -2
            switch L
                case L2-2
                    cg(n)=sqrt(((-3+L2+M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/...
                        (L2*(1-L2-4*power(L2,2)+4*power(L2,3))))/2;
                case L2-1
                    cg(n)=-(sqrt(((1+L2-M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/....
                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3))))/sqrt(2));
                case L2
                    cg(n)=sqrt(1.5)*sqrt(((1+L2-M2)*(2+L2-M2)*(-1+L2+M2)*(L2+M2))/...
                        (L2*(-3+L2*(1+4*L2*(2+L2)))));
                case L2+1
                    cg(n)=-(sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(L2+M2))/...
                        (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3))))/sqrt(2));
                case L2+2
                    cg(n)=(power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(4+L2-M2))/...
                        (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4))))/2;
            end
        case -1
            switch L
                case L2-2
                    cg(n)=-sqrt(((L2-M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/...
                        (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg(n)=((1+L2-2*M2)*sqrt(((-1+L2+M2)*(L2+M2))/...
                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2
                    cg(n)=sqrt(1.5)*sqrt(((1+L2-M2)*(L2+M2))/...
                        (L2*(-3+L2*(1+4*L2*(2+L2)))))*(-1+2*M2);
                case L2+1
                    cg(n)=-((power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2))/...
                        (L2*(1+L2)*(2+L2)*(1+2*L2)))*(L2+2*M2))/sqrt(2));
                case L2+2
                    cg(n)=power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(1+L2+M2))/...
                        (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4)));
            end
        case 0
            switch L
                case L2-2
                    cg(n)=sqrt(1.5)*sqrt(((-1+L2-M2)*(L2-M2)*(-1+L2+M2)*(L2+M2))/...
                        (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg(n)=sqrt(3)*M2*sqrt((power(L2,2)-power(M2,2))/...
                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3))));
                case L2
                    cg(n)=(-(L2*(1+L2))+3*power(M2,2))/...
                        sqrt(L2*(-3+L2*(1+4*L2*(2+L2))));
                case L2+1
                    cg(n)=-(power(-1,-2*L2+2*M2)*sqrt(3)*M2*sqrt(((1+L2-M2)*(1+L2+M2))/...
                        (L2*(1+L2)*(2+L2)*(1+2*L2))));
                case L2+2
                    cg(n)=power(-1,-2*L2+2*M2)*sqrt(3)*sqrt(((1+L2-M2)*(2+L2-M2)*(1+L2+M2)*(2+L2+M2))/...
                        (12+50*L2+70*power(L2,2)+40*power(L2,3)+8*power(L2,4)));
            end
        case 1
            switch L
                case L2-2
                    cg(n)=-sqrt(((-2+L2-M2)*(-1+L2-M2)*(L2-M2)*(L2+M2))/...
                        (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg(n)=-((sqrt(((-1+L2-M2)*(L2-M2))/...
                        ((-1+L2)*L2*(1+L2)*(1+2*L2)))*(1+L2+2*M2))/sqrt(2));
                case L2
                    cg(n)=-(power(-1,-2*L2+2*M2)*sqrt(1.5)*sqrt(((L2-M2)*(1+L2+M2))/...
                        (L2*(-3+L2*(1+4*L2*(2+L2)))))*(1+2*M2));
                case L2+1
                    cg(n)=(power(-1,-2*L2+2*M2)*(L2-2*M2)*sqrt(((1+L2+M2)*(2+L2+M2))/...
                        (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2+2
                    cg(n)=power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(1+L2+M2)*(2+L2+M2)*(3+L2+M2))/...
                        (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4)));
            end
        case 2
            switch L
                case L2-2
                    cg(n)=sqrt(((-3+L2-M2)*(-2+L2-M2)*(-1+L2-M2)*(L2-M2))/...
                        (L2*(1-L2-4*power(L2,2)+4*power(L2,3))))/2;
                case L2-1
                    cg(n)=(power(-1,-2*L2+2*M2)*sqrt(((-2+L2-M2)*(-1+L2-M2)*(L2-M2)*(1+L2+M2))/...
                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2
                    cg(n)=power(-1,-2*L2+2*M2)*sqrt(1.5)*sqrt(((-1+L2-M2)*(L2-M2)*(1+L2+M2)*(2+L2+M2))/...
                        (L2*(-3+L2*(1+4*L2*(2+L2)))));
                case L2+1
                    cg(n)=(power(-1,-2*L2+2*M2)*sqrt(((L2-M2)*(1+L2+M2)*(2+L2+M2)*(3+L2+M2))/...
                        (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2+2
                    cg(n)=(power(-1,-2*L2+2*M2)*sqrt(((1+L2+M2)*(2+L2+M2)*(3+L2+M2)*(4+L2+M2))/...
                        (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4))))/2;
            end
    end
    
end

end

% Consistency enforcement
function grumble(max_rank)
if (numel(max_rank)~=1)||(~isnumeric(max_rank))||(~isreal(max_rank))||...
   (max_rank<1)||(mod(max_rank,1)~=0)
    error('max_rank must be a real positive integer.');
end
end

% To invent, you need a good imagination and a pile of junk.
% 
% Thomas Edison

