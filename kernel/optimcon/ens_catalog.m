% Ensemble case catalog for optimal control problems. Enumerates the
% Cartesian product of the state-target pairs, the drift generators,
% the control power levels, the resonance offsets, the phase cycle
% lines, and the distortion functions; then applies the ensemble
% correlation filters and the ensemble budget. Each row of the cata-
% log is one ensemble case to be simulated at each evaluation of the
% control sequence fidelity. Syntax:
%
%              [catalog,ens_sizes,case_weights]=ens_catalog(control)
%
% Parameters:
%
%    control    - control data structure produced by optimcon.m
%
% Outputs:
%
%    catalog    - [n_cases x 6] array of ensemble indices; the col-
%                 umns index the state-target pair, the drift gene-
%                 rator, the power level, the offset combination,
%                 the phase cycle line, and the distortion function
%
%    ens_sizes  - [1 x 6] array of the ensemble dimension sizes the
%                 catalog was built from, in the same column order
%
%    case_weights - normalised [n_cases x 1] weights in catalog order
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ens_catalog.m>

function [catalog,ens_sizes,case_weights]=ens_catalog(control)

% Check consistency
grumble(control);

% Get offset ensemble size
off_ens_sizes=cellfun(@numel,control.offsets);
if ~isempty(off_ens_sizes)
    n_offset_vals=prod(off_ens_sizes);
else
    n_offset_vals=1;
end

% Extract ensemble grid dimensions
if isfield(control,'gate_targets')
    n_state_pairs=numel(control.gate_targets);
else
    n_state_pairs=numel(control.rho_init);
end
n_ens_systems=control.ndrifts;
n_power_levls=numel(control.pwr_levels);
n_phase_specs=size(control.phase_cycle,1);
n_distortions=size(control.distortion,1);

% Record the grid dimensions
ens_sizes=[n_state_pairs n_ens_systems n_power_levls ...
           n_offset_vals n_phase_specs n_distortions];

% Create a catalog of the ensemble
catalog=(1:n_state_pairs)';
catalog=[kron(ones(n_ens_systems,1),catalog) kron((1:n_ens_systems)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_power_levls,1),catalog) kron((1:n_power_levls)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_offset_vals,1),catalog) kron((1:n_offset_vals)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_phase_specs,1),catalog) kron((1:n_phase_specs)',ones(size(catalog,1),1))];
if ismember('dist_drift',control.ens_corrs)

    % Match the drift and distortion without forming unavailable pairs
    catalog=sortrows([catalog catalog(:,2)],6);

else

    % Enumerate the independent distortion dimension
    catalog=[kron(ones(n_distortions,1),catalog) kron((1:n_distortions)',ones(size(catalog,1),1))];

end

% Ensemble correlation: own state pair for each member
if ismember('rho_ens',control.ens_corrs)
    catalog=catalog(:,2:end);
    catalog=unique(catalog,'rows','stable');
    catalog=[(1:size(catalog,1))' catalog];
end

% Ensemble correlation: own state pair for each drift
if ismember('rho_drift',control.ens_corrs)
    catalog(catalog(:,1)~=catalog(:,2),:)=[];
end

% Ensemble correlation: own control power for each drift
if ismember('power_drift',control.ens_corrs)
    catalog(catalog(:,3)~=catalog(:,2),:)=[];
end

% Count the full ensemble size
n_cases=size(catalog,1);

% Convert fractional budget into sample count
ens_budget=control.budget;
if isfinite(ens_budget)&&(ens_budget<=1)
    ens_budget=round(n_cases*ens_budget);
    ens_budget=max(1,ens_budget);
end

% Apply ensemble budget
if ens_budget<n_cases

    % Stochastic selection is not a deterministic weighted quadrature
    if isfield(control,'drift_weights')
        error('Spinach:optimcon:weightedBudget',...
              'control.drift_weights requires a full ensemble budget.');
    end

    % Get RNG into a reproducible state
    rng_state=rng; rng(5318008,'twister');

    % Draw a random subset of the ensemble
    catalog=catalog(randperm(n_cases,ens_budget),:);

    % Release RNG
    rng(rng_state);

end

% Give each retained case its drift weight, then normalise once
if isfield(control,'drift_weights')

    % Refuse unequal replication without a specified joint measure
    counts=accumarray(catalog(:,2),ones(size(catalog,1),1),[control.ndrifts 1]);
    if any(counts~=counts(1))
        error('Spinach:optimcon:unequalDriftReplication',...
              'drift weights require equal retained case counts per drift.');
    end

    % Scale before summation to avoid overflow for large finite weights
    case_weights=control.drift_weights(catalog(:,2));
    scale=max(case_weights);
    if isempty(scale)||(scale==0)
        error('Spinach:optimcon:zeroDriftWeight',...
              'the retained ensemble must have positive total drift weight.');
    end
    case_weights=case_weights/scale;
    case_weights=case_weights/sum(case_weights);

else

    % Give equal weight to the retained members, including budgeted subsets
    n_cases=size(catalog,1);
    case_weights=ones(n_cases,1)/n_cases;

end

end

% Consistency enforcement
function grumble(control)
if ~isstruct(control)
    error('control must be a data structure produced by optimcon.m');
end
if ~all(isfield(control,{'ndrifts','pwr_levels','offsets',...
                         'phase_cycle','distortion','ens_corrs','budget'}))
    error('control must be a data structure produced by optimcon.m');
end
end

% The purpose of abstraction is not to be vague, but to create
% a new semantic level in which one can be absolutely precise.
%
% Edsger W. Dijkstra

