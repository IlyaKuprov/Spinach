% Strips Spinach data object down to the minimum required for a parallel
% execution stage. Syntax:
%
%             parfor_ss=stripper(spin_system,stage)
%
% Parameters:
%
%    spin_system - Spinach data object
%
%    stage       - character string specifying the parallel stage:
%
%                  'report'            - console output only
%
%                  'basis'             - basis descriptor construction:
%                                        sys.output, comp.nspins,
%                                        comp.isotopes
%
%                  'basis_projection'  - sphten-to-Zeeman projector:
%                                        bas.basis, comp.mults
%
%                  'symmetry'          - basis permutation table:
%                                        bas.basis, comp.nspins
%
%                  'operator'          - operator and state generation:
%                                        sys, tols, bas.formalism,
%                                        selected bas and comp fields
%
%                  'propagator'        - propagators and directional
%                                        derivatives: sys, tols,
%                                        bas.formalism
%
%                  'step'              - state propagation steps
%
%                  'evolution'         - independent subspace evolution:
%                                        sys, tols, bas.formalism,
%                                        selected bas fields
%
%                  'context'           - MAS rotor-stack workers
%
%                  'equilibrium'       - thermal equilibrium workers:
%                                        sys, tols, bas.formalism,
%                                        rlx.temperature, comp.mults
%
%                  'redfield_integral' - asynchronous Redfield integral
%
%                  'kinetics'          - chemical kinetics generators
%
% Outputs:
%
%    parfor_ss   - reduced Spinach data object containing only the fields
%                  required by the specified parallel stage
%
% Note: intentionally small structures are copied whole, for example
%       sys and tols; large optional structures, such as mesh and control,
%       are copied only for stages that need them.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=stripper.m>

function parfor_ss=stripper(spin_system,stage)

% Check consistency
grumble(spin_system,stage);

% Start with an empty object
parfor_ss=struct();

% Populate fields required by the specific stage
switch stage

    case 'report'

        % Keep only the output destination
        parfor_ss.sys.output=spin_system.sys.output;

    case 'basis'

        % Keep reporting and spin labels used by basis filters
        parfor_ss.sys.output=spin_system.sys.output;
        parfor_ss.comp.nspins=spin_system.comp.nspins;
        parfor_ss.comp.isotopes=spin_system.comp.isotopes;

    case 'basis_projection'

        % Keep spherical tensor basis and spin multiplicities
        parfor_ss.bas.basis=spin_system.bas.basis;
        parfor_ss.comp.mults=spin_system.comp.mults;

    case 'symmetry'

        % Keep the data needed to permute basis descriptors
        parfor_ss.bas.basis=spin_system.bas.basis;
        parfor_ss.comp.nspins=spin_system.comp.nspins;

    case {'operator','state'}

        % Keep structural data required by human2opspec and superop
        parfor_ss.sys=spin_system.sys;
        parfor_ss.tols=spin_system.tols;
        parfor_ss.bas.formalism=spin_system.bas.formalism;
        parfor_ss.comp.nspins=spin_system.comp.nspins;
        parfor_ss.comp.isotopes=spin_system.comp.isotopes;
        parfor_ss.comp.mults=spin_system.comp.mults;
        parfor_ss.comp.types=spin_system.comp.types;
        if isfield(spin_system.comp,'iso_hash')
            parfor_ss.comp.iso_hash=spin_system.comp.iso_hash;
        end

        % Keep basis tables only when they exist
        if isfield(spin_system.bas,'basis')
            parfor_ss.bas.basis=spin_system.bas.basis;
        end
        if isfield(spin_system.bas,'basis_hash')
            parfor_ss.bas.basis_hash=spin_system.bas.basis_hash;
        end
        if isfield(spin_system.bas,'lpst')
            parfor_ss.bas.lpst=spin_system.bas.lpst;
        end
        if isfield(spin_system.bas,'rpst')
            parfor_ss.bas.rpst=spin_system.bas.rpst;
        end
        if isfield(spin_system,'chem')
            parfor_ss.chem=spin_system.chem;
        end

    case {'propagator','step','evolution','redfield_integral'}

        % Keep numerical tolerances, run switches, and formalism data
        parfor_ss.sys=spin_system.sys;
        parfor_ss.tols=spin_system.tols;
        parfor_ss.bas.formalism=spin_system.bas.formalism;

        % Keep approximation selector for Hilbert-space evolution branches
        if strcmp(stage,'evolution')&&isfield(spin_system.bas,'approximation')
            parfor_ss.bas.approximation=spin_system.bas.approximation;
        end

        % Keep multiplicities for unit operators and Liouville dimensions
        if isfield(spin_system,'comp')&&isfield(spin_system.comp,'mults')
            parfor_ss.comp.mults=spin_system.comp.mults;
        end
        if isfield(spin_system.bas,'basis')
            parfor_ss.bas.basis=spin_system.bas.basis;
        end

        % Keep irrep projectors for trajectory-level symmetry reduction
        if strcmp(stage,'evolution')&&isfield(spin_system.bas,'irrep')
            parfor_ss.bas.irrep=spin_system.bas.irrep;
        end

    case 'context'

        % Keep data used by rotor-stack context workers
        parfor_ss.sys=spin_system.sys;
        parfor_ss.tols=spin_system.tols;
        parfor_ss.bas=spin_system.bas;
        parfor_ss.comp=spin_system.comp;
        parfor_ss.inter=spin_system.inter;
        if isfield(spin_system,'rlx')
            parfor_ss.rlx=spin_system.rlx;
        end
        if isfield(spin_system,'chem')
            parfor_ss.chem=spin_system.chem;
        end

    case 'equilibrium'

        % Keep thermodynamic constants and formalism dimensions
        parfor_ss.sys=spin_system.sys;
        parfor_ss.tols=spin_system.tols;
        parfor_ss.bas.formalism=spin_system.bas.formalism;
        parfor_ss.comp.mults=spin_system.comp.mults;
        parfor_ss.rlx.temperature=spin_system.rlx.temperature;

    case 'kinetics'

        % Keep reaction mapping and basis descriptors
        parfor_ss.sys.output=spin_system.sys.output;
        parfor_ss.bas.basis=spin_system.bas.basis;
        parfor_ss.chem.parts=spin_system.chem.parts;

    otherwise

        % Complain and bomb out
        error('unknown parallel execution stage.');

end

end

% Consistency enforcement
function grumble(spin_system,stage)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if (~ischar(stage))||(~ismember(stage,{'report','basis','operator',...
       'basis_projection','symmetry','state','propagator','step',...
       'evolution','context','equilibrium','redfield_integral',...
       'kinetics'}))
    error('stage must be a valid character string.');
end
end

% Insanity in individuals is something rare - but in groups, parties,
% nations, and epochs, it is the rule.
%
% Friedrich Nietzsche
