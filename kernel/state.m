% Generates Hilbert space density matrices and Liouville space state 
% vectors from their human-readable descriptions. Syntax:
%
%              rho=state(spin_system,states,spins,method)
%
% Parameters:
%
% 1. If states is a string and spins is a string
%
%                      states='Lz'; spins='13C';
%
% the function returns the sum of the corresponding single-spin densi-
% ty matrices (Hilbert space) or state vectors (Liouville space) on 
% all spins of that type. Valid labels for states in this type of call
% are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
% ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
% 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
% lid labels for spins are standard isotope names, as well as 'elect-
% rons', 'nuclei', and 'all'.
%
% 2. If states is a string and spins is a vector
%
%                      states='Lz'; spins=[1 2 4];
%
% the function returns the sum of all single-spin density matrices
% (Hilbert space) or state vectors (Liouville space) for all spins
% with the specified numbers. Valid labels for states are the same as
% in Item 1 above.
%
% 3. If states is a cell array of strings and spins is a cell array
%    of numbers:
%
%                    states={'Lz','L+'}; spins={1,2};
%
% then a product state density matrix (Hilbert space) or state vector
% (Liouville space) is produced. In the case above, Spinach will gene-
% rate LzS+ density matrix in Hilbert space or its state vector in Li-
% ouville space. Valid labels for operators are the same as in Item 1
% above.
%
% 4. For wavefunction formalism, states must be specified as an array
%    of projection quantum numbers on all spins; in that case only two
%    arguments are needed, for example, in a {'1H','1H','14N'} system:
%
%                   psi=state(spin_system,[-1/2 1/2 0])
%
% Method argument has the following effect in sphten-liouv formalism:
%
%    'cheap'  - the state vector is generated without
%               normalisation. For very large spin sys-
%               tems this is much faster
%
%    'exact'  - exact state vector with correct normalisation,
%               this is the default when the last argument is
%               skipped in the function call
%
%    'chem'   - deprecated alias for 'exact', accepted for one release
%
% Every density-matrix and Liouville method weights each substance block
% by chem.concs. Use coil_state for unweighted detection operators.
% Single-substance wavefunctions remain unweighted; segmented wavefunction
% requests are rejected (use coil_state for unweighted ket storage). The method is ignored in
% Zeeman Hilbert and Liouville formalisms, but concentration weighting is not.
%
% Outputs:
%
%     rho     - a Hilbert space density matrix or a Liouville
%               space state vector
%
% d.savostyanov@soton.ac.uk
% luke.edwards@ucl.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=state.m>

function rho=state(spin_system,states,spins,method)

% Preserve the established optional arguments
if ~exist('method','var'), method='exact'; end
if ~exist('spins','var'), spins=[]; end

% Check the formalism and wrapper-specific option
grumble(spin_system,states,spins,method);

% Retain the retired keyword for one release
if strcmp(method,'chem')
    warning('Spinach:state:deprecatedChem',...
            '''chem'' is deprecated: use state for weighted states and coil_state for unweighted coils.');
    method='exact';
end

% Construct the unweighted operator representation
rho=coil_state(spin_system,states,spins,method);

% Keep storage-only wavefunctions normalised independently of concentration
if strcmp(spin_system.bas.formalism,'zeeman-wavef'), return; end

% Weight each substance without dividing by any concentration
for n=1:spin_system.bas.nsubst
    rows=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    rho(rows,:)=spin_system.chem.concs(n)*rho(rows,:);
end

end

% Input validation function
function grumble(spin_system,states,spins,method)

if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('basis set information is missing, run basis() before calling this function.');
end
if ~ischar(spin_system.bas.formalism)
    error('formalism specification must be a character string.');
end
if ~ismember(spin_system.bas.formalism,{'zeeman-hilb', 'zeeman-liouv',...
                                        'sphten-liouv','zeeman-wavef'})
    error('unknown formalism specification.');
end
if strcmp(spin_system.bas.formalism,'zeeman-wavef')&&...
   spin_system.bas.nsubst>1
    error('Spinach:state:segmentedZeeman',...
          'concentration-weighted states are not supported in segmented zeeman-wavef formalism.');
end

if ~ischar(method)
    error('method must be a character string.')
elseif ~ismember(method, {'cheap', 'exact', 'chem'})
    error('unknown method specification.');
end

if (~(ischar(states)&&ischar(spins)))&&...
   (~(iscell(states)&&iscell(spins)))&&...
   (~(ischar(states)&&isnumeric(spins)))&&...
   (~(isnumeric(states)&&isempty(spins)))
    error('invalid state specification.');
end
if isnumeric(states)&&(numel(states)~=spin_system.comp.nspins)
    error('numel(states) must match number of spins in the system.')
end
if isnumeric(states)&&any(mod(states,0.5)~=0,'all')
    error('spin projection numbers must be integer or half-integer.');
end
if isnumeric(states)
    spin_qn=(spin_system.comp.mults(:).'-1)/2;
    if any(abs(states(:).')>spin_qn)||any(mod(states(:).'+spin_qn,1)~=0)
        error('each projection quantum number must be an allowed level of its spin.');
    end
end
if iscell(states)&&iscell(spins)&&(numel(states)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(states)&&any(~cellfun(@ischar,states))
    error('all elements of the operators cell array should be strings.');
end
if isnumeric(spins)&&(~isempty(spins))
    if (~isreal(spins))||(~isrow(spins))||any(mod(spins,1)~=0)||any(spins<1)
        error('when numeric, spins must be a row of positive integers.');
    end
    if numel(spins)~=numel(unique(spins))
        error('spin list must not have any repetitions.');
    end
end
if iscell(spins)
    if isempty(spins)
        error('when a cell array, spin list cannot be empty.');
    end
    for n=1:numel(spins)
        if (~isreal(spins{n}))||(mod(spins{n},1)~=0)||(spins{n}<1)
            error('when a cell array, spins must contain positive integers.');
        end
    end
    spins=cell2mat(spins(:));
    if numel(spins)~=numel(unique(spins))
        error('spin list must not have any repetitions.');
    end
end
end

% Aggressive public displays of virtue are where 
% the morally deplorable hide.
%
% Milo Yiannopoulos

