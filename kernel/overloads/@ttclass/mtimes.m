% Performs tensor train multiplication followed by a shrink. Syntax:
%
%                           c=mtimes(a,b)
%
% where first operand (a) can be scalar or tensor train, second operand
% (b) can be scalar, or tensor train, or full matrix.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/mtimes.m>

function c=mtimes(a,b)

% Decide the type combination
if isa(a,'ttclass')&&isa(b,'double')&&isscalar(b)
    
    % Multiply tensor train by a scalar from the right
    c=a; c.coeff=b*c.coeff; c.tolerance=b*c.tolerance;
    
elseif isa(a,'ttclass')&&isa(b,'double')&&numel(b)>1
    
    % Read sizes and ranks of the operands
    [a_ncores,a_ntrains]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a); b_sizes=size(b);

    % Check consistency
    if (prod(a_sizes(:,2))~=b_sizes(1))
        error('tensor train and vector sizes should be consistent for multiplication.');
    end
    
    % Preallocate result
    mc=prod(a_sizes(:,1)); nc=b_sizes(2); c=zeros(mc,nc);
    
    % Loop over the buffers of the operands
    for na=1:a_ntrains
        
        % Set current vector as the right-hand side
        d=b;

        % Loop over the cores and perform multiplication
        for ka=a_ncores:(-1):1
            
            % Reshape the current vector
            d=reshape(d,[a_ranks(ka+1,na)*a_sizes(ka,2),prod(a_sizes(1:ka-1,2))*prod(a_sizes(ka+1:a_ncores,1))*nc]);

            % Extract the core of the first operand and reshape it appropriately
            a_core=a.cores{ka,na};
            a_core=reshape(a_core,[a_ranks(ka,na),a_sizes(ka,1),a_sizes(ka,2),a_ranks(ka+1,na)]);
            a_core=permute(a_core,[1,2,4,3]);
            a_core=reshape(a_core,[a_ranks(ka,na)*a_sizes(ka,1),a_ranks(ka+1,na)*a_sizes(ka,2)]);

            % Perform the contraction
            d=a_core*d;
            d=reshape(d,[a_ranks(ka,na),a_sizes(ka,1),prod(a_sizes(1:ka-1,2))*prod(a_sizes(ka+1:a_ncores,1)),nc]);
            d=permute(d,[1,3,2,4]);
            
        end
        
        % Accumulate the result
        c=c+a.coeff(na)*reshape(d,[mc,nc]);
    end
       
elseif isa(a,'double')&&isa(b,'ttclass')&&isscalar(a)
    
    % Multiply tensor train by a scalar from the left
    c=b; c.coeff=a*c.coeff; c.tolerance=a*c.tolerance;
    
elseif isa(a,'ttclass')&&isa(b,'ttclass')
    
    % Read sizes and ranks of the operands
    [a_ncores,a_ntrains]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a);
    [b_ncores,b_ntrains]=size(b.cores); b_ranks=ranks(b); b_sizes=sizes(b);
    
    % Check consistency
    if (a_ncores~=b_ncores)||(~all(a_sizes(:,2)==b_sizes(:,1)))
        error('tensor train structures must be consistent.');
    end
    
    % Preallocate the result
    new_cores=cell(a_ncores,a_ntrains,b_ntrains);
    new_coeff=zeros(a_ntrains,b_ntrains);
    new_tolerance=zeros(a_ntrains,b_ntrains);
    
    % Loop over the buffers of the operands
    for nb=1:b_ntrains
        for na=1:a_ntrains
            
            % Loop over the cores
            for nc=1:a_ncores
                
                % Extract cores from the operands
                core_of_a=a.cores{nc,na};
                core_of_b=b.cores{nc,nb};
                
                % Run the index contraction: [la, ma, na, ra] * [lb, mb, nb, rb] -> [la, lb, ma, nb, ra, rb]
                core_of_a=reshape(permute(core_of_a,[1 4 2 3]),[a_ranks(nc,na)*a_ranks(nc+1,na)*a_sizes(nc,1),a_sizes(nc,2)]);
                core_of_b=reshape(permute(core_of_b,[2 3 1 4]),[b_sizes(nc,1),b_sizes(nc,2)*b_ranks(nc,nb)*b_ranks(nc+1,nb)]);
                core_of_c=reshape(core_of_a*core_of_b,[a_ranks(nc,na),a_ranks(nc+1,na),a_sizes(nc,1),b_sizes(nc,2),b_ranks(nc,nb),b_ranks(nc+1,nb)]);
                core_of_c=permute(core_of_c,[1 5 3 4 2 6]);
                
                % Store the core of the result
                new_cores{nc,na,nb}=reshape(core_of_c,[a_ranks(nc,na)*b_ranks(nc,nb),a_sizes(nc,1),b_sizes(nc,2),a_ranks(nc+1,na)*b_ranks(nc+1,nb)]);
                
            end
            
            % Multiply coefficients
            new_coeff(na,nb)=a.coeff(na)*b.coeff(nb);
            if a.coeff(na)==0 || b.coeff(nb)==0
                
                % The summand is exactly zero
                new_tolerance(na,nb)=0;
            
            else
                
                % Compute the norm of the current summand
                tt_current=ttclass;
                tt_current.coeff=a.coeff(na)*b.coeff(nb);
                tt_current.cores=cell(a_ncores,1);
                tt_current.cores=new_cores(:,na,nb);
                tt_current.tolerance=0;
                norm_current=norm(tt_current,'fro');
                
                % Update accuracy
                new_tolerance(na,nb)=norm_current*(a.tolerance(1,na)/abs(a.coeff(1,na))+b.tolerance(1,nb)/abs(b.coeff(1,nb)));
                
            end
            
        end
    end
    
    % Write the output data structure
    c=ttclass; c.coeff=reshape(new_coeff,[1,a_ntrains*b_ntrains]);
    c.cores=reshape(new_cores,[a_ncores,a_ntrains*b_ntrains]);
    c.tolerance=reshape(new_tolerance,[1,a_ntrains*b_ntrains]);
    
    % Filter out zero coeff
    pos=find(c.coeff);
    if ~isempty(pos)
        c.coeff=c.coeff(pos);
        c.cores=c.cores(:,pos);
        c.tolerance=c.tolerance(pos);
    else
        c=0*unit_like(c);
    end
    
    % Compress the answer
    c=shrink(c);
    
else
    
    % Complain and bomb out
    error('both operands must be either scalars or tensor trains.');
    
end

end

% Nothing is more revolting than the majority; for it consists of few vigorous
% predecessors, of knaves who accommodate themselves, of weak people who assi-
% milate themselves, and the mass that toddles after them  without knowing in
% the least what it wants.
%
% Johann Wolfgang von Goethe

