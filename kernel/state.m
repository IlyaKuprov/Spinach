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
% Every method weights each substance block by chem.concs. Use coil_state
% for unweighted detection operators. The method is ignored in Zeeman
% Hilbert and Liouville formalisms, but concentration weighting is not.
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

% Check the wrapper-specific option
grumble(method);

% Retain the retired keyword for one release
if strcmp(method,'chem')
    warning('Spinach:state:deprecatedChem',...
            '''chem'' is deprecated: use state for weighted states and coil_state for unweighted coils.');
    method='exact';
end

% Construct the unweighted operator representation
rho=coil_state(spin_system,states,spins,method);

% Weight each substance without dividing by any concentration
for n=1:spin_system.bas.nsubst
    rows=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    rho(rows,:)=spin_system.chem.concs(n)*rho(rows,:);
end

end

% Input validation function
function grumble(method)
if ~ischar(method)
    error('method must be a character string.');
end
end

% Aggressive public displays of virtue are where 
% the morally deplorable hide.
%
% Milo Yiannopoulos

