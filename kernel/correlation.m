% Correlation order selection function - keeps only the specified orders
% of spin correlation in the state vector. This is useful as an analyti-
% cal replacement for complicated phase cycles. Syntax:
%
%       rho=correlation(spin_system,rho,correlation_orders,spins)
%
% Parameters:
%
%   rho     -  a state vector or a horizontal stack thereof;
%              in zeeman-hilb, a density matrix or a horizon-
%              tal stack thereof
%
%   correlation_orders  -  a row vector of correlation
%                          orders to keep
%
%   spins               -  which spins to consider (e.g.
%                          '1H', '13C', 'all')
%
% Outputs:
%
%   rho     - the state vector with the undesired orders of
%             spin correlations zeroed out
%
% Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
%       hilb formalism; Fokker-Planck direct products are supported in
%       the Liouville space formalisms. In the Zeeman formalisms the
%       selection is an exact projection built from per-spin identity
%       component channels because correlation order is not diagonal
%       in the Zeeman basis; zeeman-hilb density matrices are stretch-
%       ed into Liouville space, filtered there, and folded back.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=correlation.m>

function rho=correlation(spin_system,rho,orders,spins)

% Set the default to all spins
if ~exist('spins','var'), spins='all'; end

% Check consistency
grumble(spin_system,rho,orders,spins)

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    spn_dim=spn_dim^2;
end
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Parse the spin specification
if ~isnumeric(spins)
    if strcmp(spins,'all')
        spins=1:length(spin_system.comp.isotopes);
    else
        spins=find(strcmp(spins,spin_system.comp.isotopes));
    end
end

% Decide how to proceed
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Compute the order of correlation for each basis state
        orders_present=sum(logical(spin_system.bas.basis(:,spins)),2);

        % Wipe all correlation orders except those specified by the user
        state_mask=false(size(spin_system.bas.basis,1),1);
        for n=orders
            state_mask=state_mask|(orders_present==n);
        end

        % Apply the mask
        rho(~state_mask,:)=0;

    case {'zeeman-liouv','zeeman-hilb'}

        % Hilbert space dimension and multiplicities
        dim=prod(spin_system.comp.mults);
        mults=spin_system.comp.mults;

        % Identity component channels of the selected spins
        idch=cell(1,numel(spins));
        for n=1:numel(spins)
            idch{n}=sparse(dim^2,dim^2);
            for p=1:mults(spins(n))
                for q=1:mults(spins(n))
                    A=kron(kron(speye(prod(mults(1:(spins(n)-1)))),...
                                sparse(p,q,1,mults(spins(n)),mults(spins(n)))),...
                           speye(prod(mults((spins(n)+1):end))));
                    idch{n}=idch{n}+kron(conj(A),A)/mults(spins(n));
                end
            end
        end

        % Generating function samples on the roots of unity
        nsel=numel(spins); gsam=cell(1,nsel+1);
        for j=0:nsel
            gs=rho; x=exp(2i*pi*j/(nsel+1));
            for n=1:nsel
                gs=idch{n}*gs+x*(gs-idch{n}*gs);
            end
            gsam{j+1}=gs;
        end

        % Discrete Fourier weights of the correlation orders to keep
        wts=zeros(nsel+1,1); okords=orders(ismember(orders,0:nsel));
        for j=0:nsel
            wts(j+1)=sum(exp(-2i*pi*okords(:)*j/(nsel+1)))/(nsel+1);
        end

        % Assemble the projected state
        rho=wts(1)*gsam{1};
        for j=1:nsel
            rho=rho+wts(j+1)*gsam{j+1};
        end

end

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end
    
end

% Consistency enforcement
function grumble(spin_system,rho,correlation_orders,spins)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})
    error('analytical correlation order selection is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.');
end
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if (~isnumeric(correlation_orders))||(~isvector(correlation_orders))||...
     any(correlation_orders<0)||any(mod(correlation_orders,1)~=0)
    error('correlation_orders parameter must be a vector of non-negative integers.');
end
if ~isnumeric(spins)&&~ischar(spins)
    error('spins parameter must be ''all'', an isotope label, or a list of spin numbers.');
end
if isnumeric(spins)
    if (~isreal(spins))||(~isvector(spins))||any(spins<1)||...
       any(spins>spin_system.comp.nspins)||any(mod(spins,1)~=0)
        error('spins parameter must be a vector of valid spin numbers.');
    end
elseif (~strcmp(spins,'all'))&&(~ismember(spins,spin_system.comp.isotopes))
    error('spins parameter must be ''all'' or an isotope label present in the system.');
end
end

% The first iteration of the Spin Dynamics course (http://spindynamics.org) 
% was so difficult that every single student had dropped out by about Lecture
% 10. The rest of the course was read to Rusty, Dusty, Scratchy, Patchy and
% Scruffy, the five plastic chairs in IK's office - they made for an excel-
% lent (if a bit shy so far as questions were concerned) audience. To this
% day, the total number of students who verifiably understood the whole of 
% that course can be counted on the fingers of one hand.

