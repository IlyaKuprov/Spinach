% Central-window peak-to-peak amplitude of saturated CW ESR spectra of a
% nitroxide radical as a function of the microwave field amplitude, and
% its fit to the empirical saturation law. Uses the model of Figure 2
% and compares with Table 1 of
%
%            https://doi.org/10.1007/s00723-021-01390-7
%
% The spectra are steady states of the stochastic Liouville equation
% with the microwave term retained to all orders in B1, computed with
% the magnetic parameters quoted in the caption of that figure. The
% empirical law is Equation (24) of the paper,
%
%            I(B1)=I0*B1/(1+P*B1^2)^eps
%
% for which the published fit gives P=7.263 1/G^2, eps=1.582, and a
% coefficient of determination of 0.9994. The amplitude constant I0
% is not comparable, because it carries the detection normalisation
% of the simulation rather than that of the paper. The B1 sampling is
% logarithmic and the optimiser is fminsearch, rather than the paper
% fitting protocol; P and eps are therefore independent estimates.
% The peak-to-peak estimator uses the fixed central search window; at
% strong driving its extrema may reach a window boundary as lines overlap.
%
% Conventions. Spinach refers ESR offsets to the free-electron Larmor
% frequency, and with the opposite sign to the resonance frequency: a
% g factor above the free-electron value therefore resonates at a ne-
% gative offset, and the magnetic field ascends with the offset grid.
% The field axis built below follows that convention, which places the
% broadest hyperfine line at high field, as a nitroxide requires. The
% absorption is minus the imaginary part of the steady-state transverse
% magnetisation, and the derivative is taken with respect to the field
% in Gauss, so that it is the quantity recorded by field modulation.
%
% In the isolated unsaturated Bloch limit, this model has Lorentzian
% absorption FWHM R2/pi Hz, with R1=1/T1=2*We and R2=1/T2=We.
% Additional SLE motional broadening is not a single Lorentzian width.
%
% Calculation time: seconds to minutes, depending on parallel resources
%
% alexey.bogdanov@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function sle_esr_sat_pwr()

% Static magnetic field, 3360 Gauss at X band
sys.magnet=0.3360;

% One electron and one nitrogen nucleus
sys.isotopes={'E','14N'};

% Nitroxide g tensor, diagonal in the magnetic frame
inter.zeeman.eigs={[2.0087 2.0057 2.0021],[]};
inter.zeeman.euler={[0 0 0],[]};

% Nitrogen hyperfine tensor, 6, 6, and 36 Gauss, collinear with g
inter.coupling.eigs=cell(2);
inter.coupling.euler=cell(2);
inter.coupling.eigs{1,2}=1e6*gauss2mhz([6 6 36]);
inter.coupling.euler{1,2}=[0 0 0];

% Lindblad relaxation theory in the laboratory frame
inter.relaxation={'lindblad'};
inter.rlx_keep='labframe';

% Electron rates implied by We=1e5 1/s, with no pure dephasing
inter.lind_r1_rates=[2e5 0];
inter.lind_r2_rates=[1e5 0];

% Thermalise the relaxation superoperator towards Boltzmann equilibrium
inter.equilibrium='IME';
inter.temperature=298;

% Complete spherical tensor Liouville space
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Orientational basis truncation
parameters.max_rank=10;

% Irradiated and detected isotope
parameters.spins={'E'};

% Transverse detection state
parameters.coil=coil_state(spin_system,'L+','E');

% Centre of the nitroxide spectrum on the Spinach offset axis
g_iso=mean(inter.zeeman.eigs{1});
free_g=spin_system.tols.freeg;
larmor_e=abs(spin('E'))*sys.magnet/(2*pi);
centre_freq=-(g_iso-free_g)/free_g*larmor_e;

% Frequency interval corresponding to the displayed 3300 to 3420 Gauss
parameters.sweep=centre_freq+1e6*gauss2mhz(1,g_iso)*[-60 60];
parameters.npoints=512;

% Magnetic field axis in Gauss, ascending with the offset grid
freq_axis=linspace(parameters.sweep(1),parameters.sweep(2),parameters.npoints)';
field_axis=sys.magnet*1e4+(freq_axis-centre_freq)/(1e6*gauss2mhz(1,g_iso));

% Isotropic rotational diffusion at 1e8 1/s
parameters.tau_c=1/(6*1e8);

% Microwave field amplitudes, spanning the saturation turnover
b1_grid=logspace(log10(0.02),log10(2),20);

% Preallocate the central line peak-to-peak amplitudes
ptp_amp=zeros(size(b1_grid));

% Central line search window, a tenth of the sweep either side of centre
centre_idx=round(parameters.npoints/2);
half_window=round(0.1*parameters.npoints);
region_idx=(centre_idx-half_window):(centre_idx+half_window);

% Loop over the microwave field amplitudes
for n=1:numel(b1_grid)

    % Electron nutation frequency for this B1
    parameters.irr_powers=abs(spin('E'))*b1_grid(n)*1e-4/(2*pi);

    % Steady state of the driven stochastic Liouville equation
    spectrum=gridfree(spin_system,@slowpass_sat,parameters,'esr');

    % Absorption, differentiated with respect to the magnetic field
    derivative=gradient(-imag(spectrum),field_axis);

    % Peak-to-peak amplitude of the central hyperfine line
    ptp_amp(n)=max(derivative(region_idx))-min(derivative(region_idx));

end

% Fit positive coefficients through their logarithms
sat_model=@(p,b1)exp(p(1))*b1./(1+exp(p(2))*b1.^2).^exp(p(3));

% Least squares objective over the simulated amplitudes
objective=@(p)sum((sat_model(p,b1_grid)-ptp_amp).^2);

% Refine from the published fit parameters
guess=log([max(ptp_amp) 7.263 1.582]);
p_opt=fminsearch(objective,guess,optimset('TolX',1e-10,'TolFun',1e-14,...
                                          'MaxFunEvals',1e4,'MaxIter',1e4));

% Fitted parameters and the coefficient of determination
fit_par=exp(p_opt);
residuals=sat_model(p_opt,b1_grid)-ptp_amp;
r_squared=1-sum(residuals.^2)/sum((ptp_amp-mean(ptp_amp)).^2);

% Report the fit beside the published values
disp(['saturation parameter P = ' num2str(fit_par(2))]);
disp(['correction exponent    = ' num2str(fit_par(3))]);
disp(['determination R^2      = ' num2str(r_squared)]);

% Plot the simulated amplitudes against the fitted saturation law
b1_fine=logspace(log10(min(b1_grid)),log10(max(b1_grid)),200);
kfigure('Position',[100 100 1000 760]);
plot(b1_grid,ptp_amp,'o'); hold on;
plot(b1_fine,sat_model(p_opt,b1_fine),'-');
kxlabel('B_1, Gauss');
kylabel('central-window peak-to-peak amplitude, a.u.');
klegend('SLE simulation','empirical law','Location','Best'); box on;

end


