% EPR orientation-selected quadrature grid. Syntax:
%
%   [selection_grid,selection]=epr_orisel_grid(spin_system,parameters)
%
% parameters uses the inputs of epr_orisel_transitions(), plus:
%
%   selection_transition - electron-mS row indices or 'all' (default)
%   grid_cutoff          - fraction of integrated weight to discard,
%                          in [0,1); default zero
%   plot_orisel          - plot raw selection if nonzero (optional)
%   export_grid          - MAT-file path for selected grid (optional)
%
% Outputs:
%
%   selection_grid       - alphas, betas, gammas, weights: column vectors
%                          suitable for parameters.grid in powder()
%   selection            - transition diagnostics, original grid in
%                          grid, raw response in selection_weights
%                          (also weights), selected_weights row,
%                          transition_indices, and selection_grid
%
% Raw excitation weights are multiplied by the original quadrature
% weights ONCE, sorted, truncated, and normalised to sum to one. The
% smallest prefix carrying at least 1-grid_cutoff is retained. Zero
% weight orientations are omitted. Signed quadratures are unsupported.
% An empty excitation window is an error, not a uniform powder average.
%
% The existing pruning-only call is also supported:
%
%   selection_grid=epr_orisel_grid(sph_grid,parameters)
%
% Here sph_grid contains equal-length angle/weight columns, and
% parameters.selection_weights is a row of raw excitation weights;
% parameters.grid_cutoff controls pruning. These weights must not
% already include the original quadrature weights.
%
% ilya.kuprov@weizmann.ac.il
% alexey.bogdanov@weizmann.ac.il

function [selection_grid,selection]=epr_orisel_grid(spin_system,parameters)

% Resolve established optional parameters
if ~isfield(parameters,'grid_cutoff'), parameters.grid_cutoff=0; end
if ~isfield(parameters,'selection_transition'), parameters.selection_transition='all'; end
grumble(spin_system,parameters);

% Obtain the raw selection weights
if isfield(parameters,'selection_weights')
    orisel_grid=spin_system;
    selected_weights=parameters.selection_weights;
    selection.selected_weights=selected_weights;
else
    [orisel_grid,data]=epr_orisel_transitions(spin_system,parameters);
    sel_idx=parameters.selection_transition;
    if isempty(sel_idx)||isequal(sel_idx,'all')
        sel_idx=1:size(data.selection_weights,1);
    end
    if any(sel_idx>size(data.selection_weights,1))
        error('selection transition indices are out of range.');
    end
    selected_weights=sum(data.selection_weights(sel_idx,:),1);
    selection=data;
    selection.grid=orisel_grid;
    selection.selection_weights=data.selection_weights;
    selection.weights=data.selection_weights;
    selection.selected_weights=selected_weights;
    selection.transition_indices=sel_idx;
end

% Form a non-negative quadrature measure
weights=orisel_grid.weights.*selected_weights';
if any(~isfinite(weights))||any(weights<0)
    error('orientation selection requires finite non-negative quadrature weights.');
end
if ~any(weights>0)
    error('no EPR excitation weight remains on the supplied grid.');
end
weights=weights/sum(weights);

% Retain the smallest prefix meeting the integral cutoff
[weights,idx]=sort(weights,'descend');
n=nnz(weights);
if parameters.grid_cutoff>0
    % Enforce the exact unit endpoint despite cumulative roundoff
    cumulative=cumsum(weights(1:n)); cumulative(end)=1;
    n=find(cumulative>=1-parameters.grid_cutoff,1);
end
idx=idx(1:n); weights=weights(1:n);
selection_grid.alphas=orisel_grid.alphas(idx);
selection_grid.betas=orisel_grid.betas(idx);
selection_grid.gammas=orisel_grid.gammas(idx);
selection_grid.weights=weights/sum(weights);
selection.selection_grid=selection_grid;

% Preserve optional presentation and export for the full selection call
if ~isfield(parameters,'selection_weights')
    if isfield(parameters,'plot_orisel')&&(parameters.plot_orisel~=0)
        plot_orisel(orisel_grid,selected_weights,1);
    end
    if isfield(parameters,'export_grid')&&ischar(parameters.export_grid)
        save(parameters.export_grid,'-struct','selection_grid');
    end
end

end

% Check input consistency
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isstruct(parameters)
    error('spin_system and parameters must be structures.');
end
if ~isnumeric(parameters.grid_cutoff)||~isreal(parameters.grid_cutoff)||...
   ~isscalar(parameters.grid_cutoff)||~isfinite(parameters.grid_cutoff)||...
   (parameters.grid_cutoff<0)||(parameters.grid_cutoff>=1)
    error('parameters.grid_cutoff must be a finite real scalar in [0,1).');
end
idx=parameters.selection_transition;
if ~(ischar(idx)&&strcmp(idx,'all'))&&...
   (~isnumeric(idx)||~isreal(idx)||(~isempty(idx)&&~isrow(idx))||...
    any(~isfinite(idx))||any(mod(idx,1)~=0)||any(idx<1)||...
    (numel(unique(idx))~=numel(idx)))
    error('selection_transition must contain unique positive indices or be ''all''.');
end
if isfield(parameters,'selection_weights')
    w=parameters.selection_weights;
    if ~isnumeric(w)||~isreal(w)||~isrow(w)||any(~isfinite(w))||any(w<0)
        error('parameters.selection_weights must be a finite non-negative row.');
    end
    fields={'alphas','betas','gammas','weights'};
    for k=1:numel(fields)
        x=spin_system.(fields{k});
        if ~isnumeric(x)||~isreal(x)||~iscolumn(x)||...
           (numel(x)~=numel(w))||any(~isfinite(x))
            error('grid fields must be finite real columns matching the weights.');
        end
    end
end
end


