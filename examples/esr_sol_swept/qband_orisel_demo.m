% Orientation selection of a rigid MTSSL-14N nitroxide at 34 GHz.
% Three field positions sample different parts of its Q-band spectrum.
% Colours show relative excitation intensity, normalised independently
% at each field; the returned grids have unit integrated weight.
%
% The aligned g=[2.0090 2.0060 2.0020] and A=[14 14 95] MHz tensors
% are illustrative nitroxide parameters, not a fit. This static model
% omits strain, motion, relaxation, and pulse saturation. Microwave
% excitation has a Gaussian intensity FWHM of 100 MHz, not a pulse
% amplitude bandwidth. Increase the spherical grid for convergence.
%
% Calculation time: seconds to minutes on a CPU.
%
% Written by Talos for Alexey Bogdanov

function qband_orisel_demo()

% Spin system and reference field
sys.isotopes={'E','14N'};
sys.magnet=1.21;
sys.disable={'hygiene'};
sys.parallel={'local',1};

% Aligned electron g and hyperfine tensors
inter.zeeman.eigs={[2.0090 2.0060 2.0020],[0 0 0]};
inter.zeeman.euler={[0 0 0],[0 0 0]};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([14 14 95])*1e6;
inter.temperature=20;

% Hilbert-space basis
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Microwave excitation and powder quadrature
parameters.grid='leb_2ang_rank_95';
parameters.selection_spin='E';
parameters.exclude_from_EPR=[];
parameters.selection_transition='all';
parameters.mw_frq=34e9;
parameters.bandwidth_EPR=100e6;
parameters.grid_cutoff=1e-2;
parameters.rspt_order=inf;
parameters.intensity_cutoff=0;
fields=[1.2097 1.2112 1.2162];

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Compute and display the three orientation selections
figure(); tiledlayout(1,3);
for n=1:numel(fields)
    parameters.magnetic=fields(n);
    [selection_grid,selection]=epr_orisel_grid(spin_system,parameters);
    nexttile();
    plot_orisel(selection.grid,selection.selected_weights,1);
    title(sprintf('%.4f T (%d orientations)',fields(n),numel(selection_grid.weights)));
end

end


