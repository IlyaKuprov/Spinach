% Sets up interaction tensors under partial ordering in a liquid 
% crystal with the user-supplied order matrix. All adjustable pa-
% rameters are set during the call to create.m function. Syntax:
%
%               spin_system=residual(spin_system)
%
% Parameters:
%
%      spin_system  -  the output of create.m containing
%                      spin system and interaction infor-
%                      mation, which must include the or-
%                      der matrix
%
% Outputs:
%
%      spin_system  -  the same object with anisotropic 
%                      parts of all interaction tensors
%                      replaced with their partial order
%                      residuals.
%
% Note: this function is only applicable to weak residual order
%       in high-field NMR spectroscopy.
%
% Note: the function overwrites the interaction tensors supplied
%       by the user. Relaxation superoperator, if required, must
%       be computed before this function is called.
%
% Note: this function is invoked automatically by liquid.m con-
%       text when when parameters.needs cell array contains 'rdc'.
%
% ledwards@cbs.mpg.de
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=residual.m>

function spin_system=residual(spin_system)

% Check consistency
grumble(spin_system);

% Loop over chemical subsystems
for s=1:numel(spin_system.chem.parts)

    % Process Zeeman interactions
    for n=spin_system.chem.parts{s}

        % Process non-empty interaction tensors
        if ~isempty(spin_system.inter.zeeman.matrix{n})

            % Obtain isotropic part
            iso=eye(3)*trace(spin_system.inter.zeeman.matrix{n})/3;

            % Calculate residual order
            extra_zz=trace(spin_system.inter.order_matrix{s}*(spin_system.inter.zeeman.matrix{n}-iso));

            % Update Zeeman tensor
            spin_system.inter.zeeman.matrix{n}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);

        end

    end

    % Process spin-spin couplings
    pairs=spin_system.inter.coupling.matrix;
    indices=find(ismember(cell2mat(pairs(:,1)),spin_system.chem.parts{s}));
    drop=false(size(pairs,1),1);
    for q=indices'

        % Average the anisotropy over the residual order tensor
        iso=trace(pairs{q,3})*eye(3)/3;
        extra_zz=trace(spin_system.inter.order_matrix{s}*(pairs{q,3}-iso));
        pairs{q,3}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);
        drop(q)=norm(pairs{q,3},2)<2*pi*spin_system.tols.inter_cutoff;

    end
    spin_system.inter.coupling.matrix=pairs(~drop,:);
    if isfield(spin_system.inter.coupling,'strength')
        spin_system.inter.coupling.strength(drop,:)=[];
    end

end

% Report back to the user
report(spin_system,'interaction tensors have been replaced by their weak order residuals.');

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system.inter,'order_matrix'))||isempty(spin_system.inter.order_matrix)
    error('order matrix infomation is missing from the spin_system structure.');
end
end

% If I'd observed all the rules, I'd never have got anywhere.
%
% Marilyn Monroe


