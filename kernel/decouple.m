% Obliterates all interactions and populations in the subspace of states
% that involve the specified spins in any way. The specified spins would
% not contribute to the system dynamics until the Liouvillian is rebuilt
% from scratch. Syntax:
%
%                 [L,rho]=decouple(spin_system,L,rho,spins)
%
% Parameters:
%
%     L     - Liouvillian superoperator, this may be left empty
%
%     rho   - state vector or a horizontal stack thereof, this
%             may be left empty
%
%     spins - spins to be wiped, specified either by name, e.g.
%             {'13C','1H'}, or by a list of numbers, e.g. [1 2]
%
% Outputs:
%
%     rho   - state vector(s) with all populations of the 
%             states involving the target spins set to zero
%
%     L     - Liouvillian superoperator with all terms in-
%             volving the target spins set to zero 
%
% Note: this function is an analytical equivalent of a perfect decoup-
%       ling pulse sequence on the specified spins.
%
% Note: this function requires sphten-liouv formalism and supports Fok-
%       ker-Planck direct products.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=decouple.m>

function [L,rho]=decouple(spin_system,L,rho,spins)

% Return if the spin list is empty
if isempty(spins), return; end

% Check consistency
grumble(spin_system,L,rho,spins);

% Find the spins to be decoupled
if isnumeric(spins)
    
    % Find spins by numbers
    dec_mask=false(1,spin_system.comp.nspins); dec_mask(spins)=true;
    
else
    
    % Find spins by name
    dec_mask=ismember(spin_system.comp.isotopes,spins);
    
end

% Inform the user
report(spin_system,[num2str(nnz(dec_mask)) ' spins to be frozen and depopulated.']);
                
% Get the list of states to be wiped
zero_mask=(sum(spin_system.bas.basis(:,dec_mask),2)~=0);

% Process the Liouvillian
if (nargout>0)&&(~isempty(L))
    
    % Get dimension statistics
    spn_dim=size(spin_system.bas.basis,1); 
    spc_dim=size(L,1)/spn_dim;
    
    % Kron the list into the Fokker-Planck space
    fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask));

    % Inform the user
    report(spin_system,['space sub-problem dimension: ' num2str(spc_dim)]);
    report(spin_system,['spin sub-problem dimension:  ' num2str(spn_dim)]);
    report(spin_system,['zeroing ' num2str(nnz(fp_zero_mask))...
                        ' rows and columns in the Liouvillian.']);

    % Apply the zero mask
    L(fp_zero_mask,:)=0; L(:,fp_zero_mask)=0;
    
    % Re-evaluate sparsity
    L=clean_up(spin_system,L,spin_system.tols.liouv_zero);

end

% Process state vector stack
if (nargout>1)&&(~isempty(rho))
    
    % Get dimension statistics
    spn_dim=size(spin_system.bas.basis,1); 
    spc_dim=size(rho,1)/spn_dim;

    % Kron the list into the Fokker-Planck space
    fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask));
    
    % Inform the user
    report(spin_system,['space sub-problem dimension: ' num2str(spc_dim)]);
    report(spin_system,['spin sub-problem dimension:  ' num2str(spn_dim)]);
    report(spin_system,['zeroing ' num2str(nnz(fp_zero_mask))...
                        ' rows in the state vector.']);
                    
    % Apply the zero mask
    rho(fp_zero_mask,:)=0;
    
end

end

% Consistency enforcement
function grumble(spin_system,L,rho,spins)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('analytical decoupling is only available for sphten-liouv formalism.');
end
if (~isempty(rho))&&(~isempty(L))&&(size(L,2)~=size(rho,1))
    error('matrix dimensions of L and rho must agree.');
end
if size(L,1)~=size(L,2)
    error('Liouvillian must be a square matrix.');
end
if (~isnumeric(spins))&&(~iscell(spins))
    error('spins parameter must either be a list of numbers or a cell array of strings.');
end
if iscell(spins)&&any(~ismember(spins,spin_system.comp.isotopes))
    error('the system does not contain the spins specified.');
end
if isnumeric(spins)&&(size(spins,1)~=1)
    error('if spins are specified by number, a row vector of numbers must be used.');
end
if isnumeric(spins)&&(max(spins)>spin_system.comp.nspins)
    error('the spin number specified is greater than the number of spins in the system.');
end
if isnumeric(spins)&&(any(~isreal(spins))||any(spins<1))
    error('spin numbers must be real positive integers.');
end
end

% It's not worth doing something unless you were doing something that
% someone, somewere, would much rather you weren't doing.
%
% Terry Pratchett 

