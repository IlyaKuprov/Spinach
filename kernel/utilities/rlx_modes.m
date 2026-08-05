% Bosonic mode dissipation superoperator. Builds thermalised GKSL
% dissipators for the amplitude damping and the pure dephasing of
% the bosonic modes declared in inter.modes, using the amplitude
% damping rates and the pure dephasing rates ingested by create.m
% and the Bose-Einstein thermal occupation numbers computed from
% the physical mode frequencies, meaning the sum of the declared
% carrier and the declared frequency where inter.modes.carriers
% is present, and the system temperature. Syntax:
%
%                     R=rlx_modes(spin_system)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%                   with bosonic mode information present
%
% Outputs:
%
%    R            - bosonic mode dissipation superoperator
%
% Note: the dissipators are kappa*(1+nbar)*D[a], kappa*nbar*D[c],
%       and 2*gamma_phi*D[n], where D[x] is the GKSL dissipator of
%       the operator x, built from ladder operators truncated to
%       the level count of each mode. The Spinach convention of
%       zero temperature meaning the high-temperature limit is not
%       applicable to bosonic modes: zero temperature here produces
%       zero thermal occupation numbers.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_modes.m>

function R=rlx_modes(spin_system)

% Check consistency
grumble(spin_system);

% Preallocate the answer
R=mprealloc(spin_system,1);

% Locate bosonic modes
mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}));

% Loop over bosonic modes
for k=mode_list

    % Get the dissipation rates
    kappa=spin_system.inter.modes.damp(k);
    gphi=spin_system.inter.modes.dephase(k);

    % Move on if the mode is not dissipative
    if (kappa==0)&&(gphi==0), continue; end

    % Get the thermal occupation number, only damping needs it
    if (kappa>0)&&(spin_system.rlx.temperature>0)

        % Physical frequency, laboratory carrier included where declared
        if spin_system.inter.modes.carriers(k)>0
            phys_frq=spin_system.inter.modes.carriers(k)+...
                     spin_system.inter.modes.frqs(k);
        else
            phys_frq=abs(spin_system.inter.modes.frqs(k));
        end

        % Catch modes whose thermal occupation is undefined
        if phys_frq<2*pi*spin_system.tols.inter_cutoff
            error(['thermal occupation of zero-frequency mode ' num2str(k) ...
                   ' is undefined, supply the laboratory frame frequency.']);
        end

        % Bose-Einstein statistics at the physical frequency
        beta_factor=spin_system.tols.hbar*phys_frq/...
                    (spin_system.tols.kbol*spin_system.rlx.temperature);
        nbar=1/(exp(beta_factor)-1);

        % Inform the user which frequency has been used
        report(spin_system,['bosonic mode ' num2str(k) ': thermal occupation from a '...
                            'physical frequency of ' num2str(phys_frq/(2*pi)) ' Hz, nbar = ' num2str(nbar)]);

    else

        % Zero occupation without damping or at zero temperature
        nbar=0;

    end

    % Inform the user
    report(spin_system,['bosonic mode ' num2str(k) ': kappa = ' num2str(kappa) ...
                        ' s^-1, gamma_phi = ' num2str(gphi) ' s^-1, nbar = ' num2str(nbar)]);

    % Get the ladder superoperators
    a_left=operator(spin_system,{'A'},{k},'left');
    a_right=operator(spin_system,{'A'},{k},'right');
    c_left=operator(spin_system,{'C'},{k},'left');
    c_right=operator(spin_system,{'C'},{k},'right');
    n_left=operator(spin_system,{'N'},{k},'left');
    n_right=operator(spin_system,{'N'},{k},'right');

    % Add the cooling dissipator
    R=R+kappa*(1+nbar)*(a_left*c_right-0.5*(n_left+n_right));

    % Add the heating dissipator
    if nbar>0
        ac_left=operator(spin_system,{'AC'},{k},'left');
        ac_right=operator(spin_system,{'AC'},{k},'right');
        R=R+kappa*nbar*(c_left*a_right-0.5*(ac_left+ac_right));
    end

    % Add the pure dephasing dissipator
    if gphi>0
        R=R+2*gphi*(n_left*n_right-0.5*(n_left*n_left+n_right*n_right));
    end

end

end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system.inter,'modes')
    error('bosonic mode information is missing from the spin_system structure.');
end
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('bosonic mode dissipators are only available in Liouville space.');
end
if ~isfield(spin_system.rlx,'temperature')
    error('temperature information is missing from the spin_system structure.');
end
end

% To achieve great things, two things are needed: a plan,
% and not quite enough time.
%
% Leonard Bernstein


