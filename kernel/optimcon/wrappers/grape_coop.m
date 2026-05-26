% Pairs of cooperative pulses that may be used as components of a phase
% cycle. The pulses are designed to produce as much of the destination
% state as they can, and to have imputities of opposite sign. Adding the
% outcomes of the two experiments then destroys the impurities. Syntax:
%
%   [traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)
%
% Parameters:
%
%      phi_profile  -  phase profiles of the two pulses,
%                      concatenated horizontally
%
% Outputs:
%
%      traj_data    -  trajectory information structure
%
%      fidelity     -  cooperative fidelity measure
%
%      gradient     -  cooperative fidelity gradient
%
% Note: only phase-modulated point-to-point transformations are supported.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=grape_coop.m>

function [traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)

% Check consistency
grumble(spin_system);

% Extract phase profiles
n_channels=spin_system.control.ncontrols/2;
profile_a=phi_profile(1:n_channels,:);
profile_b=phi_profile((n_channels+1):end,:);

% Get target and impurity projectors
rho_targ=spin_system.control.rho_targ{1};
switch spin_system.bas.formalism
    case 'zeeman-hilb'
        targ_norm=hdot(rho_targ,rho_targ);
        if abs(targ_norm)==0
            error('target state has zero norm.');
        end
    otherwise
        P_targ=rho_targ*rho_targ';
        P_dirt=eye(size(P_targ))-P_targ;
end

% Make sure final states are available
spin_system.control.return_traj=true();

% Preserve the caller's ensemble correlations
ens_corrs=spin_system.control.ens_corrs;

% Run both experiments
[traj_data_a,fidelity_a,gradient_a]=grape_phase(profile_a,spin_system);
[traj_data_b,fidelity_b,gradient_b]=grape_phase(profile_b,spin_system);

% Project out the impurities
dirt_a=cell(numel(traj_data_a),1);
for n=1:numel(traj_data_a)
    switch spin_system.bas.formalism
        case 'zeeman-hilb'
            rho_a=traj_data_a{n}.forward{end};
            dirt_a{n}=rho_a-rho_targ*hdot(rho_targ,rho_a)/targ_norm;
        otherwise
            dirt_a{n}=P_dirt*traj_data_a{n}.forward(:,end);
    end
end
dirt_b=cell(numel(traj_data_b),1);
for n=1:numel(traj_data_b)
    switch spin_system.bas.formalism
        case 'zeeman-hilb'
            rho_b=traj_data_b{n}.forward{end};
            dirt_b{n}=rho_b-rho_targ*hdot(rho_targ,rho_b)/targ_norm;
        otherwise
            dirt_b{n}=P_dirt*traj_data_b{n}.forward(:,end);
    end
end
dirt_sum=cell(size(dirt_a));
for n=1:numel(dirt_sum)
    dirt_sum{n}=dirt_a{n}+dirt_b{n};
end

% Get offset ensemble size
off_ens_sizes=cellfun(@numel,spin_system.control.offsets);
if ~isempty(off_ens_sizes)
    n_offset_vals=prod(off_ens_sizes);
else
    n_offset_vals=1;
end

% Extract non-rho ensemble grid dimensions
n_ens_systems=spin_system.control.ndrifts;
n_power_levls=numel(spin_system.control.pwr_levels);
n_phase_specs=size(spin_system.control.phase_cycle,1);
n_distortions=size(spin_system.control.distortion,1);

% Recreate the original non-rho ensemble catalogue
catalog=(1:n_ens_systems)';
catalog=[kron(ones(n_power_levls,1),catalog) kron((1:n_power_levls)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_offset_vals,1),catalog) kron((1:n_offset_vals)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_phase_specs,1),catalog) kron((1:n_phase_specs)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_distortions,1),catalog) kron((1:n_distortions)',ones(size(catalog,1),1))];

% Preserve own-state-per-drift filtering from the main pass
if ismember('rho_drift',ens_corrs)
    catalog(catalog(:,1)~=1,:)=[];
end

% Preserve own-power-per-drift filtering from the main pass
if ismember('power_drift',ens_corrs)
    catalog(catalog(:,2)~=catalog(:,1),:)=[];
end

% Count the full impurity ensemble size
n_cases=size(catalog,1);

% Get ensemble budget
if isfield(spin_system.control,'budget')
    ens_budget=spin_system.control.budget;
else
    ens_budget=Inf;
end

% Convert fractional budget into sample count
if isfinite(ens_budget)&&(ens_budget<=1)
    ens_budget=round(n_cases*ens_budget);
    ens_budget=max(1,ens_budget);
end

% Apply ensemble budget
if ens_budget<n_cases

    % Get RNG into a reproducible state
    rng_state=rng; rng(5318008,'twister');

    % Draw the same random subset of the ensemble
    catalog=catalog(randperm(n_cases,ens_budget),:);

    % Release RNG
    rng(rng_state);

end

% Check that each impurity target has one ensemble case
if size(catalog,1)~=numel(dirt_sum)
    error('impurity target count does not match the cooperative ensemble catalogue.');
end
spin_system.control.rho_targ=dirt_sum;

% Replicate the initial state
rho=spin_system.control.rho_init{1};
spin_system.control.rho_init=cell(size(spin_system.control.rho_targ));
spin_system.control.rho_init(:)={rho};

% Impurity cancellation gradients
spin_system.control.ens_corrs={};
spin_system.control.ens_catalog=[(1:size(catalog,1))' catalog];
spin_system.control.penalties={'none'};
spin_system.control.p_weights=0;
[~,~,gradient_c]=grape_phase(profile_a,spin_system);
[~,~,gradient_d]=grape_phase(profile_b,spin_system);

% Average fidelity of the two pulses
fidelity=(fidelity_a+fidelity_b)/2;

% Penalty on the squared norm of the dirt
switch spin_system.bas.formalism
    case 'zeeman-hilb'
        fidelity(1)=fidelity(1)-mean(cellfun(@(x)real(hdot(x,x)),dirt_sum));
    otherwise
        fidelity(1)=fidelity(1)-mean(cellfun(@(x)norm(x,2)^2,dirt_sum));
end

% Assemble the physical-fidelity gradient
gradient=cat(1,gradient_a,gradient_b)/2;
gradient(:,:,1)=gradient(:,:,1)-2*cat(1,gradient_c(:,:,1),gradient_d(:,:,1));

% Return both trajectories
traj_data={traj_data_a,traj_data_b};

end

% Consistency enforcement
function grumble(spin_system)
if (numel(spin_system.control.rho_targ)~=1)||...
   (numel(spin_system.control.rho_init)~=1)
    error('this function only supports point-to-point transformations.');
end
if mod(spin_system.control.ncontrols,2)~=0
    error('grape_coop is phase-modulated, number of controls must be even.');
end
end

% Morally authoritarian movements are attractive to
% ugly, miserable, talentless people.
%
% Milo Yiannopoulos
