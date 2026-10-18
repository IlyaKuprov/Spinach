% Point-to-point magnetisation transfer in compiled local bases. Syntax:
%
%                       F=rlx_flux(spin_system)
%
% Parameters:
%
%    spin_system - Spinach data structure with a compiled basis and
%                  validated chem.flux_rate and chem.flux_type fields
%
% Outputs:
%
%    F           - sparse magnetisation flux relaxation superoperator
%
% i.kuprov@soton.ac.uk

function F=rlx_flux(spin_system)

% Check consistency
grumble(spin_system);

% Allocate the generator and retain the no-flux default
F=mprealloc(spin_system,0);
if nnz(spin_system.chem.flux_rate)==0, return; end

% The original flux model uses spherical tensor basis labels
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('magnetisation flux treatment is only available for sphten-liouv formalism.');
end

% Select the physical fate of correlations, as in the original flux model
switch spin_system.chem.flux_type
    case 'intramolecular'
        transfer_correlations=true;
    case 'intermolecular'
        transfer_correlations=false;
    otherwise
        error('unrecognized type of magnetisation flux.');
end

% Process every matrix entry, including diagonal bath rates
[sources,destins,rates]=find(spin_system.chem.flux_rate);
orders=cellfun(@(b)sum(logical(b),2),spin_system.bas.basis,'UniformOutput',false);
for n=1:numel(rates)

    % Locate single- and multi-spin orders on either side
    spins=[sources(n) destins(n)];
    rows=cell(2,2); labels=cell(2,2);
    for side=1:2
        subst=which_subst(spin_system,spins(side));
        local=find(spin_system.chem.parts{subst}==spins(side));
        descriptor=spin_system.bas.basis{subst};
        counts=orders{subst};
        for order=1:2
            if order==1
                mask=(counts==1)&(descriptor(:,local)~=0);
            else
                mask=(counts>1)&(descriptor(:,local)~=0);
            end
            rows{side,order}=spin_system.bas.offsets(subst)+find(mask);
            labels{side,order}=descriptor(mask,local);
        end
    end

    % Single-spin orders transfer in both models; correlations only within molecules
    for order=1:2
        if (order==1)||transfer_correlations

            % Orders containing both spins were stationary in the original model
            [~,from_static,to_static]=intersect(rows{1,order},rows{2,order});
            rows{1,order}(from_static)=[]; labels{1,order}(from_static)=[];
            rows{2,order}(to_static)=[]; labels{2,order}(to_static)=[];

            % Retain the original tensor-label compatibility check
            if ~isequal(labels{1,order},labels{2,order})
                error('spin systems on either side of the reaction arrow have different topology or basis sets.');
            end

            % Gain terms retain their original ordered basis-state matching
            F=F+rates(n)*sparse(rows{2,order},rows{1,order},...
                               ones(size(rows{1,order})),size(F,1),size(F,2));
        end

        % Intermolecular flux destroys every source correlation, even on the diagonal
        F=F-rates(n)*sparse(rows{1,order},rows{1,order},...
                           ones(size(rows{1,order})),size(F,1),size(F,2));
    end
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

