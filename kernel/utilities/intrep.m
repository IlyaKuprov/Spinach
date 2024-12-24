% Interaction representation transformation with respect to 
% a specified Hamiltonian to specified order in perturbation
% theory (https://doi.org/10.1063/1.4928978). Syntax:
%
%             Hr=intrep(spin_system,H0,H,T,order)
%
% Parameters:
%
%    H0     - the Hamiltonian with respect to which the
%             interaction representation transformation
%             is to be done, typically Zeeman Hamiltonian
%
%    H      - laboratory frame Hamiltonian H0+H1 that is
%             to be transformed into the interaction rep-
%             resentation, typically the full Hamiltonian
%
%    T      - period of the H0 propagator
%
%    order  - perturbation theory order in the rotating
%             frame transformation, this may be inf
%
% Outputs:
%
%    Hr     - Hamiltonian in the interaction representation
%
% Note: the auxiliary matrix method is massively faster than
%       either commutator series or diagonalisation.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=intrep.m>

function Hr=intrep(spin_system,H0,H,T,order)

% Check consistency
grumble(H0,H,order);

% Confirm that T is indeed the period
P=propagator(spin_system,H0,T);
if norm(P-speye(size(P)),1)>1e-6
    error('the value of T supplied is not a period of exp(-i*H0*t) propagator.');
end

% Compute the norms
norm_h0=norm(H0,1); norm_h1=norm(H-H0,1);

% Report the norms
report(spin_system,['rotating frame period: ' num2str(T) ' seconds']);
report(spin_system,['H0 (lab frame) 1-norm: ' num2str(norm_h0)]);
report(spin_system,['H1 (lab frame) 1-norm: ' num2str(norm_h1)]);
    
% Confirm that norms are sensible
if norm_h1>norm_h0
    error('norm of H1 found to be bigger than norm of H0: H1 must be a perturbation.');
end

% Decide the theory order
switch order
    
    case 0
        
        % Shortcut for high field
        Hr=H-H0;
        
    case Inf
        
        % Shortcut for infinite order
        Hr=(1i/T)*logm(expm(full(-1i*H*T)));
        
    otherwise
        
        % Get the derivatives
        D=dirdiff(spin_system,H0,H-H0,T,order+1);
        
        % Get the first term
        Hr=(1i/T)*(D{1}'*D{2});
        
        % Get the rest of the series
        for n=2:order
            for k=1:n
                Hr=Hr+(1i/T)*nchoosek(n-1,k-1)*D{n-k+1}'*D{k+1}/factorial(n);
            end
        end

end

% Clean up the output
Hr=clean_up(spin_system,(Hr+Hr')/2,spin_system.tols.liouv_zero);

% Return matrix density statistics
report(spin_system,['int. rep. Hamiltonian dimension ' num2str(size(Hr,1)) ', nnz ' ...
                    num2str(nnz(Hr)) ', density ' num2str(100*nnz(Hr)/numel(Hr))...
                    '%, sparsity ' num2str(issparse(Hr))]);

end

% Consistency enforcement
function grumble(H0,H,order)
if (~ishermitian(H))||(~ishermitian(H0))
    error('both H and H0 must be Hermitian.');
end
if ((~isreal(order))||(order<1)||(mod(order,1)~=0))&&(~isinf(order))
    error('unsupported rotating frame transformation theory order.');
end
end

% There is no difference between communism and socialism, except in the
% means of achieving the same ultimate end: communism proposes to enslave
% men by force, socialism - by vote. It is merely the difference between
% murder and suicide.
%
% Ayn Rand

