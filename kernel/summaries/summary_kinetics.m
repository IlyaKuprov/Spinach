% Chemical kinetics and point-to-point magnetisation flux reporting. Syntax:
%
%                  summary_kinetics(spin_system)
%                  summary_kinetics(spin_system,reactions)
%
% Parameters:
%
%    spin_system - Spinach data structure
%
%    reactions   - optional compiled reaction records; when supplied,
%                  report the reaction network instead of flux rates
%
% Outputs:
%
%    This function prints a summary to the configured output stream.
%
% i.kuprov@soton.ac.uk

function summary_kinetics(spin_system,reactions)

% Check consistency
if nargin==1
    grumble(spin_system,{});
else
    grumble(spin_system,reactions);
end

% Report flux inputs before basis construction
if nargin==1
    if isfield(spin_system.chem,'flux_rate')&&nnz(spin_system.chem.flux_rate)>0
        switch spin_system.chem.flux_type
            case 'intramolecular'
                report(spin_system,'intramolecular magnetisation flux, correlations will be kept.');
            case 'intermolecular'
                report(spin_system,'intermolecular magnetisation flux, correlations will be damped.');
        end
        [rows,cols,vals]=find(spin_system.chem.flux_rate);
        report(spin_system,'point-to-point flux rates:');
        report(spin_system,'===============================');
        report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
        report(spin_system,'-------------------------------');
        for n=1:numel(vals)
            report(spin_system,sprintf(' %-3d      %-3d      %+0.3e',rows(n),cols(n),vals(n)));
        end
    end
else

    % Report the compiled reaction network and atom transport
    summary_chemistry(spin_system);
    for n=1:numel(reactions)
        reaction=reactions{n};
        source_spins=cellfun(@(x)x(:)',spin_system.chem.parts(reaction.reactants),'UniformOutput',false);
        source_spins=[source_spins{:}];
        report(spin_system,['reaction ' num2str(n) ': matched spins ' ...
               mat2str(reaction.matching(:,1)') ', traced spins ' ...
               mat2str(setdiff(source_spins,reaction.matching(:,1))) ...
               ', closure ' reaction.closure]);
    end
end

end

% Consistency enforcement
function grumble(spin_system,reactions)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if ~iscell(reactions)||~all(cellfun(@isstruct,reactions))
    error('reactions must be a cell array of structures.');
end
end

