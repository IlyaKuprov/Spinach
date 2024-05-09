% Computes dipolar couplings in the presence or absence of periodic
% boundary conditions using an optimized approach with the conmat function.
% This function is part of the Spinach library's kernel utilities and is
% designed to efficiently handle large systems by reducing computational
% complexity from quadratic to approximately N log N.
%
% Usage of this function within the Spinach library is internal, and
% direct calls are generally discouraged in favor of higher-level functions.
% This function updates the spin_system object with new interaction matrices
% reflecting computed dipolar interactions.
%
% Contributors:
%    Ilya Kuprov - initial implementation
%    Alessandro Lodi - optimisation using conmat (lodialessandro0@gmail.com)
%
% Contact:
%    lodialessandro0@gmail.com
%
% References:
%    TODO

function spin_system = dipolar_linear(spin_system)
    % Check consistency
    grumble(spin_system);

    % Report to the user
    report(spin_system, ['dipolar interaction distance threshold: ' ...
                        num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);
    report(spin_system, 'running linear dipolar interaction network analysis with PBC handling...');

    % Aggregate all coordinates from the cell array
    all_coords = vertcat(spin_system.inter.coordinates{:});

    % Ensure no empty cells in coordinates
    if any(cellfun(@isempty, spin_system.inter.coordinates))
        error('One or more coordinate cells are empty, ensure all are filled.');
    end

    % Handle PBC logic
    pbc_vectors = spin_system.inter.pbc;
    dd_ncells = spin_system.tols.dd_ncells;

    % Loop over chemical species
    for m=1:numel(spin_system.chem.parts)
        
        % Extract the spin list
        spin_list = spin_system.chem.parts{m};
        coords = all_coords(ismember(1:length(all_coords), spin_list), :);

        % Calculate extended coordinates for PBCs if defined
        if ~isempty(pbc_vectors)
            coords = apply_pbc(coords, pbc_vectors, dd_ncells);
        end

        % Calculate the connectivity matrix using conmat
        conmatrix = conmat(coords, spin_system.tols.prox_cutoff);

        % Process interactions based on the connectivity matrix
        [rows, cols] = find(conmatrix);
        for idx = 1:length(rows)
            n = rows(idx);
            k = cols(idx);
            dv = coords(k,:) - coords(n,:);
            distance = norm(dv);

            if distance < spin_system.tols.prox_cutoff
                % Compute the orientation vector
                ort = dv / distance;

                % Compute the dipolar interaction prefactor
                A = 0.5 * spin_system.inter.gammas(spin_list(n)) * spin_system.inter.gammas(spin_list(k)) * spin_system.tols.hbar * spin_system.tols.mu0 / (4 * pi * (distance * 1e-10)^3);

                % Compute the dipolar coupling matrix
                % D = -A * (3*(ort'*ort) - eye(3));
                
                % Compute the dipolar coupling matrix
                D = A * [1-3*ort(1)^2   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
                -3*ort(2)*ort(1)  1-3*ort(2)^2   -3*ort(2)*ort(3);
                -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)^2];

                % Scale by chemical shifts and g-tensors
                D = spin_system.inter.zeeman.ddscal{spin_list(n)}' * D * spin_system.inter.zeeman.ddscal{spin_list(k)};

                % Store the computed matrix
                if isempty(spin_system.inter.coupling.matrix{spin_list(n), spin_list(k)})
                    spin_system.inter.coupling.matrix{spin_list(n), spin_list(k)} = D;
                else
                    spin_system.inter.coupling.matrix{spin_list(n), spin_list(k)} = spin_system.inter.coupling.matrix{spin_list(n), spin_list(k)} + D;
                end
            end
        end
    end
end

function extended_coords = apply_pbc(coords, pbc, dd_ncells)
    % Apply periodic boundary conditions to extend coordinates
    % coords: original coordinates array
    % pbc: periodic boundary conditions vector(s)
    % dd_ncells: number of cells in each direction to consider
    
    % Check dimensions and prepare variables
    num_atoms = size(coords, 1);
    dim_pbc = numel(pbc);  % Dimensionality of the PBC
    pbc_range = -dd_ncells:dd_ncells;

    % Determine the number of vectors based on PBC dimensionality
    if isempty(pbc)
        extended_coords = coords; % No PBCs, return original
    elseif isscalar(pbc)
        nvecs = (2 * dd_ncells + 1);
        extended_coords = zeros(num_atoms * nvecs, 3);
    elseif dim_pbc == 2
        nvecs = (2 * dd_ncells + 1)^2;
        extended_coords = zeros(num_atoms * nvecs, 3);
    elseif dim_pbc == 3
        nvecs = (2 * dd_ncells + 1)^3;
        extended_coords = zeros(num_atoms * nvecs, 3);
    else
        error('PBC translation vector array has invalid dimensions.');
    end
    
    % Generate extended coordinates based on PBCs
    linear_index = 1;
    for i = 1:num_atoms
        if isempty(pbc)
            extended_coords(linear_index, :) = coords(i, :);
            linear_index = linear_index + 1;
        else
            for p = pbc_range
                if isscalar(pbc)
                    % 1D PBC
                    extended_coords(linear_index, :) = coords(i, :) + p * pbc{1};
                    linear_index = linear_index + 1;
                elseif dim_pbc == 2
                    for q = pbc_range
                        % 2D PBC
                        extended_coords(linear_index, :) = coords(i, :) + p * pbc{1} + q * pbc{2};
                        linear_index = linear_index + 1;
                    end
                elseif dim_pbc == 3
                    for q = pbc_range
                        for r = pbc_range
                            % 3D PBC
                            extended_coords(linear_index, :) = coords(i, :) + p * pbc{1} + q * pbc{2} + r * pbc{3};
                            linear_index = linear_index + 1;
                        end
                    end
                end
            end
        end
    end
end

% Consistency enforcement
function grumble(spin_system)
    if ~all(isfield(spin_system,{'comp','inter','chem','tols'}))
        error('spin_system object is missing essential information.');
    end
end
    


