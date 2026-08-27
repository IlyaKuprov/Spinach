% Reduced effective-model illustration of the trade-off discussed by
% Nehra, Agarwal, and Nishiyama for 14N decoupling under 1H detection.
% The example is intentionally qualitative rather than quantitative:
% it shows why low-power CW decoupling is narrowband, why increasing
% the 14N RF field produces a Bloch-Siegert shift on the observed 1H
% resonance, and why B1 inhomogeneity turns that shift into broadening.
%
% The "offset-tolerant low-power" trace represents the design goal of
% Bloch-Siegert-aware robust optimal control or low-power amplitude-
% modulated decoupling. Replace the effective coefficients below with a
% more detailed Hamiltonian model for quantitative work.
%
% aditya.dev@weizmann.ac.il

function yusuke_14n_broadening_demo()

% Magnetic field corresponding to 800 MHz 1H
magnet=18.8;
nu_n14_hz=abs(spin('14N')*magnet/(2*pi));

% Three inequivalent 14N sites around the decoupler carrier
site_offsets_hz=[-18e3 0 18e3];
site_weights=[0.30 0.40 0.30];

% B1 distribution across the sample volume
b1_scales=linspace(0.85,1.15,81);
b1_weights=gaussian_weights(b1_scales,1.0,0.06);

% Frequency axes for the plots
offset_axis_hz=linspace(-40e3,40e3,801);
freq_axis_hz=linspace(-800,800,4001);

% Effective decoupling settings
rf_low_hz=8e3;
rf_high_hz=20e3;
rf_ot_hz=12e3;

% Reduced-model parameters controlling the observed 1H line
intrinsic_fwhm_hz=80;
residual_penalty_hz=220;
bs_coeff_hz=8e-7;

% Decoupling profiles across the 14N offset range
eta_low_cw=decoupling_profile(offset_axis_hz,rf_low_hz,1.0,2);
eta_high_cw=decoupling_profile(offset_axis_hz,rf_high_hz,1.0,2);
eta_low_ot=decoupling_profile(offset_axis_hz,rf_ot_hz,2.2,4);

% Predicted 1H line shapes
[line_low,stats_low]=proton_line(freq_axis_hz,site_offsets_hz,...
    site_weights,b1_scales,b1_weights,rf_low_hz,1.0,2,...
    intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz);

[line_high_hom,stats_high_hom]=proton_line(freq_axis_hz,site_offsets_hz,...
    site_weights,1.0,1.0,rf_high_hz,1.0,2,...
    intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz);

[line_high_inhom,stats_high_inhom]=proton_line(freq_axis_hz,site_offsets_hz,...
    site_weights,b1_scales,b1_weights,rf_high_hz,1.0,2,...
    intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz);

[line_low_ot,stats_low_ot]=proton_line(freq_axis_hz,site_offsets_hz,...
    site_weights,b1_scales,b1_weights,rf_ot_hz,2.2,4,...
    intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz);

% Plot the reduced-model story
figure;

subplot(2,1,1);
plot(offset_axis_hz/1e3,eta_low_cw,'LineWidth',1.5); hold on;
plot(offset_axis_hz/1e3,eta_high_cw,'LineWidth',1.5);
plot(offset_axis_hz/1e3,eta_low_ot,'LineWidth',1.5);
for n=1:numel(site_offsets_hz)
    xline(site_offsets_hz(n)/1e3,'k:','LineWidth',0.8);
end
grid on;
xlabel('^{14}N offset / kHz');
ylabel('Effective decoupling efficiency');
title(['Reduced model at B_0 = ' num2str(magnet) ' T, \nu_{14N} = ' ...
       num2str(nu_n14_hz/1e6,'%.1f') ' MHz']);
legend({'Low-power CW','High-power CW','Low-power offset-tolerant'},...
       'Location','SouthWest');

subplot(2,1,2);
plot(freq_axis_hz,normalise_line(line_low),'LineWidth',1.5); hold on;
plot(freq_axis_hz,normalise_line(line_high_hom),'LineWidth',1.5);
plot(freq_axis_hz,normalise_line(line_high_inhom),'LineWidth',1.5);
plot(freq_axis_hz,normalise_line(line_low_ot),'LineWidth',1.5);
grid on;
xlabel('Observed ^{1}H frequency relative to mean / Hz');
ylabel('Normalised intensity');
title('Bloch-Siegert shift becomes broadening once B_1 is inhomogeneous');
legend({'Low-power CW','High-power CW, homogeneous B_1',...
        'High-power CW, inhomogeneous B_1',...
        'Low-power offset-tolerant'},'Location','NorthWest');

% Console summary
fprintf('\n');
fprintf('Yusuke reduced-model summary at %.1f T\n',magnet);
fprintf('  Low-power CW:            mean decoupling = %.3f, FWHM = %.1f Hz\n',...
        stats_low.mean_efficiency,stats_low.fwhm_hz);
fprintf('  High-power CW, hom B1:   mean decoupling = %.3f, FWHM = %.1f Hz\n',...
        stats_high_hom.mean_efficiency,stats_high_hom.fwhm_hz);
fprintf('  High-power CW, inh B1:   mean decoupling = %.3f, FWHM = %.1f Hz\n',...
        stats_high_inhom.mean_efficiency,stats_high_inhom.fwhm_hz);
fprintf('  Low-power offset-tolerant: mean decoupling = %.3f, FWHM = %.1f Hz\n',...
        stats_low_ot.mean_efficiency,stats_low_ot.fwhm_hz);
fprintf('\n');
fprintf(['Interpretation: high power broadens the observed line only once the ' ...
         'Bloch-Siegert shift is dispersed by B1 inhomogeneity.\n']);
fprintf(['The low-power offset-tolerant trace indicates the design target for ' ...
         'BSS-aware robust optimal control.\n']);

end

function eta=decoupling_profile(offset_hz,rf_hz,bandwidth_gain,profile_order)

% Simple offset-response model: wider bandwidth at larger RF field or
% more sophisticated modulation, flatter passband at higher order.
bw_hz=bandwidth_gain*rf_hz;
eta=1./(1+(abs(offset_hz)./bw_hz).^profile_order);

end

function [line_shape,stats]=proton_line(freq_axis_hz,site_offsets_hz,...
    site_weights,b1_scales,b1_weights,rf_hz,bandwidth_gain,profile_order,...
    intrinsic_fwhm_hz,residual_penalty_hz,bs_coeff_hz)

% Normalise the ensemble weights
site_weights=site_weights/sum(site_weights);
b1_weights=b1_weights/sum(b1_weights);

% Weighted mean Bloch-Siegert shift, removed to focus on broadening
mean_shift_hz=0;
for n=1:numel(site_offsets_hz)
    for k=1:numel(b1_scales)
        member_weight=site_weights(n)*b1_weights(k);
        member_shift=bs_coeff_hz*(rf_hz*b1_scales(k))^2;
        mean_shift_hz=mean_shift_hz+member_weight*member_shift;
    end
end

% Assemble the proton line
line_shape=zeros(size(freq_axis_hz));
mean_efficiency=0;
for n=1:numel(site_offsets_hz)

    % Residual broadening when the current 14N site is not well decoupled
    eta=decoupling_profile(site_offsets_hz(n),rf_hz,bandwidth_gain,profile_order);
    member_fwhm_hz=intrinsic_fwhm_hz+residual_penalty_hz*(1-eta);

    % Average decoupling efficiency over the site ensemble
    mean_efficiency=mean_efficiency+site_weights(n)*eta;

    % B1 spread turns the Bloch-Siegert shift into broadening
    for k=1:numel(b1_scales)
        member_weight=site_weights(n)*b1_weights(k);
        member_shift=bs_coeff_hz*(rf_hz*b1_scales(k))^2-mean_shift_hz;
        line_shape=line_shape+member_weight*lorentzian(freq_axis_hz,...
                   member_shift,member_fwhm_hz);
    end

end

% Summary statistics
stats.mean_efficiency=mean_efficiency;
stats.mean_shift_hz=mean_shift_hz;
stats.fwhm_hz=linewidth(freq_axis_hz,line_shape);

end

function y=lorentzian(x,centre_hz,fwhm_hz)

% Lorentzian line with unit integral
gamma=fwhm_hz/2;
y=(gamma/pi)./((x-centre_hz).^2+gamma^2);

end

function weights=gaussian_weights(grid,centre,width)

% Truncated Gaussian weights for the B1 distribution
weights=exp(-((grid-centre)/width).^2/2);
weights=weights/sum(weights);

end

function line_out=normalise_line(line_in)

% Normalise for plotting
line_out=line_in/max(line_in);

end

function fwhm_hz=linewidth(freq_axis_hz,line_shape)

% Numerical full width at half maximum
line_shape=line_shape/max(line_shape);
idx=find(line_shape>=0.5);
if isempty(idx)
    fwhm_hz=NaN;
else
    fwhm_hz=freq_axis_hz(idx(end))-freq_axis_hz(idx(1));
end

end
