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
%                      stacked in two row blocks
%
% Outputs:
%
%      traj_data    -  trajectory information for both pulses
%
%      fidelity     -  cooperative fidelity measure
%
%      gradient     -  cooperative fidelity gradient
%
% Note: only phase-modulated point-to-point transformations are supported.
% The control.freeze mask uses the same two row blocks as phi_profile;
% each pulse freezes its own phase coordinates. An empty mask freezes none.
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

% Split the input-coordinate freeze mask between the two pulses
freeze_a=spin_system.control.freeze; freeze_b=freeze_a;
if ~isempty(freeze_a)
    freeze_a=freeze_a(1:n_channels,:);
    freeze_b=freeze_b((n_channels+1):end,:);
end

% Get target and impurity projectors
rho_targ=spin_system.control.rho_targ{1};
switch spin_system.bas.formalism
    case 'zeeman-hilb'
        targ_norm=hdot(rho_targ,rho_targ);
    otherwise
        targ_norm=rho_targ'*rho_targ;
        P_dirt=eye(numel(rho_targ))-rho_targ*rho_targ'/targ_norm;
end

% Make sure final states are available
spin_system.control.return_traj=true();

% Run both experiments
spin_system.control.freeze=freeze_a;
[traj_data_a,fidelity_a,gradient_a]=grape_phase(profile_a,spin_system);
spin_system.control.freeze=freeze_b;
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
spin_system.control.rho_targ=dirt_sum;

% Replicate the initial state
rho=spin_system.control.rho_init{1};
spin_system.control.rho_init=cell(size(spin_system.control.rho_targ));
spin_system.control.rho_init(:)={rho};

% The squared impurity norm requires a real linear overlap
spin_system.control.fidelity='real';

% Impurity cancellation gradients
spin_system.control.ens_corrs={'rho_ens'};
[spin_system.control.catalog,...
 spin_system.control.ens_sizes]=ens_catalog(spin_system.control);
spin_system.control.freeze=freeze_a;
[~,~,gradient_c]=grape_phase(profile_a,spin_system);
spin_system.control.freeze=freeze_b;
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

% Assemble the gradient, impurity term only enters the fidelity slice
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
if norm(spin_system.control.rho_targ{1},'fro')==0
    error('target state has zero norm.');
end
if mod(spin_system.control.ncontrols,2)~=0
    error('grape_coop is phase-modulated, number of controls must be even.');
end
if strcmp(spin_system.control.fid_type,'average')||(~isempty(spin_system.control.traj_pen))
    error('trajectory cost terms are not available in grape_coop.');
end
end

% Morally authoritarian movements are attractive to
% ugly, miserable, talentless people.
%
% Milo Yiannopoulos


