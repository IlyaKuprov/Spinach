% Saturated CW ESR spectra of a nitroxide radical as a function of the
% rotational diffusion rate, computed as the steady state of the sto-
% chastic Liouville equation with the microwave term retained to all
% orders in B1. Illustrates the model used for Figure 3 of
%
%            https://doi.org/10.1007/s00723-021-01390-7
%
% The magnetic parameters are those quoted in the caption of that fi-
% gure: gxx=2.0087, gyy=2.0057, gzz=2.0021, Axx=Ayy=6 G, Azz=36 G, an
% electron longitudinal relaxation rate We=1e5 1/s, a static field of
% 3360 Gauss, and a microwave field amplitude of 0.5 Gauss. The four
% diffusion rates are read from the figure panel labels. Vertical
% units follow Spinach thermal magnetisation at 298 K; no independent
% peak normalisation is applied.
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
% The no-pure-dephasing Lindblad model is an explicit assumption;
% the paper does not specify its transverse-rate inputs completely.
% These curves are not asserted to be an exact numerical reproduction.
%
% In the isolated unsaturated Bloch limit, this model has Lorentzian
% absorption FWHM R2/pi Hz, with R1=1/T1=2*We and R2=1/T2=We.
% Additional SLE motional broadening is not a single Lorentzian width.
%
% Calculation time: seconds to minutes, depending on parallel resources
%
% alexey.bogdanov@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function sle_esr_sat_rate()

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
parameters.mode='saturated';
parameters.sweep=centre_freq+1e6*gauss2mhz(1,g_iso)*[-60 60];
parameters.npoints=512;

% Magnetic field axis in Gauss, ascending with the offset grid
freq_axis=linspace(parameters.sweep(1),parameters.sweep(2),parameters.npoints)';
field_axis=sys.magnet*1e4+(freq_axis-centre_freq)/(1e6*gauss2mhz(1,g_iso));

% Microwave field amplitude of 0.5 Gauss, as a nutation frequency
parameters.irr_powers=abs(spin('E'))*0.5*1e-4/(2*pi);

% Isotropic rotational diffusion rates for the four panels, 1/s
r_grid=10.^[7 7.5 8 9];

% Preallocate the derivative spectra
spectra=zeros(parameters.npoints,numel(r_grid));

% Loop over the rotational diffusion rates
for n=1:numel(r_grid)

    % Correlation time corresponding to this diffusion rate
    parameters.tau_c=1/(6*r_grid(n));

    % Steady state of the driven stochastic Liouville equation
    spectrum=gridfree(spin_system,@slowpass,parameters,'esr');

    % Absorption, differentiated with respect to the magnetic field
    spectra(:,n)=gradient(-imag(spectrum),field_axis);

end

% Plot one panel per rotational diffusion rate
kfigure('Position',[100 100 1000 760]);
for n=1:numel(r_grid)

    % Derivative spectrum on the ascending field axis
    subplot(2,2,n); plot(field_axis,spectra(:,n));
    kxlabel('magnetic field, Gauss'); kylabel('first derivative, a.u.');
    ktitle(['$R$ = ' num2str(r_grid(n),'%.4g') ' s$^{-1}$']);
    xlim([min(field_axis) max(field_axis)]); box on;

end

end


