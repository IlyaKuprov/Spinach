% EPR frequencies, intensities, and orientation-selection weights. Syntax:
%
%   [sph_grid,transitions]=epr_orisel_transitions(spin_system,parameters)
%
% Parameters:
%
%   grid               - name of a Spinach spherical grid
%   selection_spin     - one electron isotope or label, default spins{1}
%                        if supplied, otherwise 'E'
%   exclude_from_EPR   - spin indices, logical mask, isotope/label, or a
%                        cell array of these; removed before selection
%                        (default: none); selection_spin must remain
%   magnetic           - positive static field in tesla
%   mw_frq             - microwave frequency in Hz
%   bandwidth_EPR      - Gaussian FWHM in Hz of the INTENSITY response;
%                        the response has unit peak, not unit area, and
%                        is applied once to transition intensities
%   rspt_order         - 1, 2, 3, 4, or Inf (default: exact diagonalisation)
%   intensity_cutoff   - discard intensities at or below this value
%                        (default: zero)
%
% Outputs:
%
%   sph_grid           - original grid; angle and weight column vectors
%   transitions        - frequencies (Hz), intensities, transition_rows,
%                        dfdb_hz_per_t (Hz/T): cells of columns, one
%                        per orientation; selection_weights is a
%                        (2S)-by-norient array; orientation_weights
%                        is the original quadrature column; parameters
%                        contains the resolved parameters
%
% Intensities are |<i|Sx|j>|^2 times the thermal population difference.
% Rows approximate adjacent electron-mS manifolds using rounded <Sz>.
% This lab-Z assignment is only meaningful in the high-field regime;
% strongly mixed transitions can be mislabelled or omitted even when
% all rows are requested. Multiple-electron selection is unsupported.
% Selection uses one electron spin and fixed laboratory-x polarization.
% It is not a pulse simulation: saturation, relaxation and strain are
% absent. The input spin system must have a nonzero reference field.
%
% For the existing selector-only entry point, set spin_selector and
% call indices=epr_orisel_transitions(spin_system,parameters). No EPR
% Hamiltonian is built; selector_name is accepted but not used.
%
% ilya.kuprov@weizmann.ac.il
% alexey.bogdanov@weizmann.ac.il

function [sph_grid,transitions]=epr_orisel_transitions(spin_system,parameters)

% Preserve the established selector-only entry point
if isfield(parameters,'spin_selector')
    grumble(spin_system,parameters);
    sph_grid=spin_selection_indices(spin_system,parameters.spin_selector);
    transitions=[]; return
end

% Resolve established optional parameters
if ~isfield(parameters,'rspt_order'), parameters.rspt_order=inf; end
if ~isfield(parameters,'intensity_cutoff'), parameters.intensity_cutoff=0; end
if ~isfield(parameters,'exclude_from_EPR'), parameters.exclude_from_EPR=[]; end
if ~isfield(parameters,'selection_spin')
    if isfield(parameters,'spins')&&iscell(parameters.spins)&&~isempty(parameters.spins)
        parameters.selection_spin=parameters.spins{1};
    else
        parameters.selection_spin='E';
    end
end
grumble(spin_system,parameters);
parameters.exclude_from_EPR=spin_selection_indices(spin_system,parameters.exclude_from_EPR);

% Remove spins excluded from the EPR Hamiltonian
loc_spin=spin_system;
if ~isempty(parameters.exclude_from_EPR)
    loc_spin=kill_spin(loc_spin,parameters.exclude_from_EPR);
end
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
loc_spin=basis(loc_spin,bas);
selection_idx=spin_selection_indices(loc_spin,parameters.selection_spin);
if (numel(selection_idx)~=1)||~strncmp(loc_spin.comp.isotopes{selection_idx},'E',1)
    error('selection_spin must identify exactly one retained electron spin.');
end

% Build Hamiltonians at the input reference field
[i_c,q_c]=hamiltonian(assume(loc_spin,'labframe','couplings'));
[i_z,q_z]=hamiltonian(assume(loc_spin,'labframe','zeeman'));
sz=operator(loc_spin,'Lz',selection_idx);
h_mw=(operator(loc_spin,'L+',selection_idx)+...
      operator(loc_spin,'L-',selection_idx))/2;
n_rows=loc_spin.comp.mults(selection_idx)-1;

% Load the original quadrature grid
sph_grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',parameters.grid),...
              'alphas','betas','gammas','weights');
alphas=sph_grid.alphas; betas=sph_grid.betas; gammas=sph_grid.gammas;
n_orients=numel(sph_grid.weights);
selection_weights=zeros(n_rows,n_orients);
freq_cell=cell(1,n_orients); int_cell=cell(1,n_orients);
row_cell=cell(1,n_orients); dfdb_cell=cell(1,n_orients);

% Use the caller's existing parallel pool, if any
pool=gcp('nocreate'); nworkers=0;
if ~isempty(pool), nworkers=min(pool.NumWorkers,n_orients); end
parfor(n=1:n_orients,nworkers)

    % Assemble the orientation-dependent reference-field Hamiltonian
    angles=[alphas(n) betas(n) gammas(n)];
    h_z=i_z+orientation(q_z,angles); h_z=(h_z+h_z')/2;
    h_c=i_c+orientation(q_c,angles); h_c=(h_c+h_c')/2;
    [e,v,de,t,lp]=rspt_eig(loc_spin,parameters,h_z,h_c,h_mw,parameters.magnetic);

    % Enumerate positive-intensity transitions
    ne=numel(e); ntran=0;
    freq_hz=zeros(ne*(ne-1)/2,1); intensity=freq_hz;
    row_idx=freq_hz; dfdb_hz=freq_hz;
    selection_row=zeros(n_rows,1);
    sz_diag=real(diag(v'*sz*v));
    for k=1:ne-1
        for m=k+1:ne
            iij=abs(t(k,m))^2*(lp(k)-lp(m));
            if iij>parameters.intensity_cutoff
                idx=round(min(sz_diag([k m]))+n_rows/2+1);
                if (idx>=1)&&(idx<=n_rows)
                    ntran=ntran+1;
                    freq_hz(ntran)=(e(m)-e(k))/(2*pi);
                    intensity(ntran)=iij;
                    row_idx(ntran)=idx;
                    dfdb_hz(ntran)=(de(m)-de(k))/(2*pi);
                end
            end
        end
    end

    % Apply a unit-peak Gaussian intensity response
    freq_hz=freq_hz(1:ntran); intensity=intensity(1:ntran);
    row_idx=row_idx(1:ntran); dfdb_hz=dfdb_hz(1:ntran);
    response=gaussfun(freq_hz-parameters.mw_frq,parameters.bandwidth_EPR)/...
             gaussfun(0,parameters.bandwidth_EPR);
    for k=1:ntran
        selection_row(row_idx(k))=selection_row(row_idx(k))+intensity(k)*response(k);
    end
    selection_weights(:,n)=selection_row;
    freq_cell{n}=freq_hz; int_cell{n}=intensity;
    row_cell{n}=row_idx; dfdb_cell{n}=dfdb_hz;

end

% Return raw transition data and unnormalised excitation weights
transitions.frequencies=freq_cell; transitions.intensities=int_cell;
transitions.transition_rows=row_cell; transitions.dfdb_hz_per_t=dfdb_cell;
transitions.selection_weights=selection_weights;
transitions.orientation_weights=sph_grid.weights;
transitions.parameters=parameters;

end

% Resolve spin indices, masks, isotope names, labels, and their cell arrays
function indices=spin_selection_indices(spin_system,selector)

% Resolve the selector according to its type
if isempty(selector)
    indices=[];
elseif isnumeric(selector)
    indices=unique(selector(:).','stable');
elseif islogical(selector)
    indices=find(selector(:).');
elseif ischar(selector)
    indices=find(strcmp(spin_system.comp.isotopes,selector));
    if isfield(spin_system.comp,'labels')
        indices=union(indices,find(strcmp(spin_system.comp.labels,selector)),'stable');
    end
    if isempty(indices), error('spin selector did not match any spins.'); end
elseif iscell(selector)
    indices=cell(size(selector));
    for n=1:numel(selector)
        indices{n}=spin_selection_indices(spin_system,selector{n});
    end
    indices=unique([indices{:}],'stable');
else
    error('unsupported spin selector type.');
end

% Check the resolved indices and mask length
if islogical(selector)&&(numel(selector)~=spin_system.comp.nspins)
    error('logical spin selector must have one entry per spin.');
end
if ~isreal(indices)||any(~isfinite(indices))||any(mod(indices,1)~=0)||...
   any(indices<1)||any(indices>spin_system.comp.nspins)
    error('spin selector indices are out of range.');
end

end

% Check input consistency
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isstruct(parameters)
    error('spin_system and parameters must be structures.');
end
if isfield(parameters,'spin_selector'), return; end
if ~ischar(parameters.grid)||~isrow(parameters.grid)
    error('parameters.grid must be a spherical grid name.');
end
if ~ischar(parameters.selection_spin)||~isrow(parameters.selection_spin)
    error('parameters.selection_spin must be an electron isotope or label.');
end
if ~isnumeric(parameters.magnetic)||~isreal(parameters.magnetic)||...
   ~isscalar(parameters.magnetic)||~isfinite(parameters.magnetic)||...
   (parameters.magnetic<=0)||(spin_system.inter.magnet<=0)
    error('selection and reference magnetic fields must be positive.');
end
if ~isnumeric(parameters.mw_frq)||~isreal(parameters.mw_frq)||...
   ~isscalar(parameters.mw_frq)||~isfinite(parameters.mw_frq)||(parameters.mw_frq<=0)
    error('parameters.mw_frq must be a positive finite frequency in Hz.');
end
if ~isnumeric(parameters.bandwidth_EPR)||~isreal(parameters.bandwidth_EPR)||...
   ~isscalar(parameters.bandwidth_EPR)||~isfinite(parameters.bandwidth_EPR)||...
   (parameters.bandwidth_EPR<=0)
    error('parameters.bandwidth_EPR must be a positive finite FWHM in Hz.');
end
if ~isnumeric(parameters.intensity_cutoff)||~isreal(parameters.intensity_cutoff)||...
   ~isscalar(parameters.intensity_cutoff)||~isfinite(parameters.intensity_cutoff)||...
   (parameters.intensity_cutoff<0)
    error('parameters.intensity_cutoff must be a finite non-negative scalar.');
end
if ~isnumeric(parameters.rspt_order)||~isreal(parameters.rspt_order)||...
   ~isscalar(parameters.rspt_order)||~ismember(parameters.rspt_order,[1 2 3 4 inf])
    error('parameters.rspt_order must be 1, 2, 3, 4, or Inf.');
end
end


