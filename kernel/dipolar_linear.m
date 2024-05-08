function spin_system = dipolar_linear(spin_system)

    % Check consistency
    grumble(spin_system);
    
    % Report to the user
    report(spin_system,['dipolar interaction distance threshold: ' ...
                        num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);
    report(spin_system,'running dipolar interaction network analysis using conmat...');
    
    % Collect all coordinates
    all_coords = vertcat(spin_system.inter.coordinates{:});
    
    % Calculate the connectivity matrix using conmat
    conmatrix = conmat(all_coords, spin_system.tols.prox_cutoff);
    
    % Find interacting atom pairs through connectivity matrix
    [rows, cols] = find(conmatrix);
    
    % Report to the user
    report(spin_system,['found ' num2str(numel(rows)/2) ' spin pair(s) with distance under the '...
                        'threshold of ' num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);
    
    % Loop over interacting atom pairs
    for idx = 1:numel(rows)
        n = rows(idx);
        k = cols(idx);
    
        % Get the distance and orientation
        distance = norm(all_coords(n,:) - all_coords(k,:), 2);
        ort = (all_coords(k,:) - all_coords(n,:)) / distance;
    
        % Compute the dipolar interaction prefactor (0.5 due to double counting)
        A = 0.5 * spin_system.inter.gammas(n) * spin_system.inter.gammas(k) * ...
            spin_system.tols.hbar * spin_system.tols.mu0 / (4 * pi * (distance * 1e-10)^3);
    
        % Compute the dipolar coupling matrix
        D = A * [1 - 3*ort(1)*ort(1),   -3*ort(1)*ort(2),   -3*ort(1)*ort(3);
                 -3*ort(2)*ort(1),  1 - 3*ort(2)*ort(2),   -3*ort(2)*ort(3);
                 -3*ort(3)*ort(1),   -3*ort(3)*ort(2),  1 - 3*ort(3)*ort(3)];
    
        % Account for chemical shifts and g-tensors
        D = spin_system.inter.zeeman.ddscal{n}' * D * spin_system.inter.zeeman.ddscal{k};
    
        % Clean up numerical noise
        D = D - eye(3) * trace(D) / 3;
    
        % Add to the total
        if isempty(spin_system.inter.coupling.matrix{n, k})
            spin_system.inter.coupling.matrix{n, k} = D;
        else
            spin_system.inter.coupling.matrix{n, k} = ...
            spin_system.inter.coupling.matrix{n, k} + D;
        end
    end
    
    end
    