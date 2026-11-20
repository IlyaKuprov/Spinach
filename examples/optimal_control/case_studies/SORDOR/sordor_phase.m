% Phase-only SORDOR design at a fixed quadratic dispersion setting.
% Scores the unmodified BURBOP seed before a short optimisation. This is
% a reduced offset grid and not the published dispersion-parameter sweep.
% Source: https://doi.org/10.1016/j.jmro.2026.100218
%
% Talos <talos@spindynamics.org>

function sordor_phase()

% Prepare a single proton with two transverse control axes
sys.isotopes={'1H'};
sys.magnet=1;
sys.parallel={'processes',2};
sys.output='hush';
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=basis(create(sys,struct()),bas);
S=pauli(2);

% Read the source waveform without changing its amplitudes or phases
seed=readmatrix(fullfile(fileparts(mfilename('fullpath')),'burbop_seed.txt'));
control.amplitudes=seed(:,1).'/1e4;
phase=deg2rad(seed(:,2).');
control.pulse_dt=0.5e-6*ones(size(phase));
control.pwr_levels=2*pi*1e4;
control.isotopes={'1H'};
control.channels=[1;1];
control.operators={S.x,S.y};
control.method='lbfgs';
control.max_iter=2;
control.penalties={'none'};
control.p_weights=0;

% Match each drift offset to its frequency-dependent quarter turn
offsets=linspace(-20e3,20e3,9);
dispersion=pi*(sum(control.pulse_dt)/25e-6)*0.25*(1-(offsets/40e3).^2);
control.ens_corrs={'rho_drift'};
control.drifts=cell(size(offsets));
control.gate_targets=cell(size(offsets));
for n=1:numel(offsets)

    % Construct and check the complete target rotation
    control.drifts{n}={2*pi*offsets(n)*S.z};
    axis=[cos(dispersion(n));sin(dispersion(n));0];
    target=expm(-1i*(pi/2)*(axis(1)*S.x+axis(2)*S.y));
    controls={S.x,S.y,S.z};
    rotation=zeros(3);
    for k=1:3
        for j=1:3
            rotation(j,k)=real(hdot(controls{j},target*controls{k}*target'))/...
                          real(hdot(controls{j},controls{j}));
        end
    end
    cross_axis=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
    assert(norm(rotation-axis*axis'-cross_axis,'fro')<1e-12);
    control.gate_targets{n}=target;

end
assert(isequal(dispersion,fliplr(dispersion)));

% Score the unchanged seed and refine its phases
spin_system=optimcon(spin_system,control);
[~,before]=grape_phase(phase,spin_system);
phase=fmaxnewton(spin_system,@grape_phase,phase);
[~,after]=grape_phase(phase,spin_system);
fprintf('SORDOR seed and refined gate fidelity: %.10g %.10g\n',before(1),after(1));

end


