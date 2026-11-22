% Computes dipolar couplings in the presence or absence of periodic
% boundary conditions. This is an auxiliary function of Spinach ker-
% nel, direct calls are discouraged. Use xyz2dd and xyz2hfc to con-
% vert Cartesian coordinates into dipolar and hyperfine couplings 
% respectively. Syntax:
%
%                  spin_system=dipolar(spin_system)
%
% Parameters:
%
%    spin_system   - Spinach data object containing infor-
%                    mation about chemical subsystems, ato-
%                    mic coordinates, and periodic bounda-
%                    ry conditions
%
% Outputs:
%
%    spin_system   - Spinach data object with the interac-
%                    tion arrays updated with dipolar and
%                    hyperfine coupling information
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dipolar.m>

function spin_system=dipolar(spin_system)

% Check consistency
grumble(spin_system);

% Report to the user
report(spin_system,['dipolar interaction distance threshold: ' ...
                    num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);
report(spin_system,'running dipolar interaction network analysis...');

% Enumerate translations in the order used for the tensor sum
axes=spin_system.inter.pbc;
shifts=zeros(1,0);
for d=1:numel(axes)
    shells=(-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells)';
    shifts=[repelem(shifts,numel(shells),1) repmat(shells,size(shifts,1),1)];
end

% Collect displacement vectors within each chemical species
blocks=cell(numel(spin_system.chem.parts),size(shifts,1));
for m=1:numel(spin_system.chem.parts)

    % Omit atoms without coordinates
    spins=spin_system.chem.parts{m};
    spins=spins(~cellfun(@isempty,spin_system.inter.coordinates(spins)));
    xyz=vertcat(spin_system.inter.coordinates{spins});
    if isempty(xyz)||(spin_system.tols.prox_cutoff==0), continue; end

    % Search the reference cell against each translated cell
    for t=1:size(shifts,1)
        images=xyz;
        for d=1:numel(axes)
            images=images+shifts(t,d)*axes{d};
        end
        [from,to]=nearby_points(xyz,images,spin_system.tols.prox_cutoff);
        dv=images(to,:)-xyz(from,:);
        keep=(from~=to)&(vecnorm(dv,2,2)<spin_system.tols.prox_cutoff);
        keys=[from(keep) to(keep)];
        blocks{m,t}=[reshape(spins(keys),size(keys)) dv(keep,:)];
    end

end

% Group all retained images of each ordered pair
hits=vertcat(blocks{:});
prox=cell(0,3);
if ~isempty(hits)
    hits=sortrows(hits,[2 1]);
    starts=[1;find(any(diff(hits(:,1:2),1,1),2))+1];
    lengths=diff([starts;size(hits,1)+1]);
    prox=[num2cell(hits(starts,1:2)) mat2cell(hits(:,3:5),lengths,3)];
end
spin_system.inter.proximity=prox;
rows=cell2mat(prox(:,1));
cols=cell2mat(prox(:,2));
tensors=prox;

% Report to the user
report(spin_system,['found ' num2str(numel(rows)/2) ' spin pair(s) with distance under the '...
                    'threshold of ' num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);

% Loop over interacting atom pairs
for n=1:numel(rows)
    
    % Accumulate couplings from individual periodic images
    tensors{n,3}=zeros(3);

    % Loop over PBC directions
    for k=1:size(prox{n,3},1)
        
        % Get the distance
        distance=norm(prox{n,3}(k,:),2);
        
        % Reject atomic collisions
        if distance<0.5
            error('collision detected between spin %d and spin %d or their PBC images.',rows(n),cols(n));
        end

        % Get the ort
        ort=prox{n,3}(k,:)/distance;
        
        % Compute the dipolar interaction prefactor (0.5 due to double counting)
        A=0.5*spin_system.inter.gammas(rows(n))*spin_system.inter.gammas(cols(n))*...
              spin_system.tols.hbar*spin_system.tols.mu0/(4*pi*(distance*1e-10)^3);
        
        % Compute the dipolar coupling matrix
        D=A*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
              -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
              -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
          
        % Optional approximate spin-orbit correction
        if ismember('sodd',spin_system.sys.enable)
            D=spin_system.inter.zeeman.ddscal{rows(n)}'*D*...
              spin_system.inter.zeeman.ddscal{cols(n)};
        end
      
        % Clean up numerical noise
        D=D-eye(3)*trace(D)/3;
        
        % Add to the total
        tensors{n,3}=tensors{n,3}+D;

    end
    
end

% Combine point dipoles with any existing tensors
spin_system.inter.coupling.matrix=sum_couplings([spin_system.inter.coupling.matrix;tensors]);

end

% Consistency enforcement
function grumble(spin_system)
if ~all(isfield(spin_system,{'comp','inter','chem','tols'}))
    error('spin_system object is missing essential information.');
end
end

% What worries me about religion is that it teaches people to
% be satisfied with not understanding. 
%
% Richard Dawkins


