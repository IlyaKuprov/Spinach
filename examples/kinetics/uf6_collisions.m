% Collision-induced quadrupolar relaxation in a simplified UF6 model.
% Each chemical species contains one representative 19F-235U pair,
% not the full six-fluorine molecule. A symmetric state (Cq=0) exchanges
% with a collision-distorted state (Cq=400 MHz); J=213 Hz in both.
% The forward rate varies while the distorted-state lifetime stays at
% 1 ps. Stationary chemical populations weight the initial coherence.
%
% Isotropic Redfield relaxation uses tau_c=10 fs at every rate. This
% phenomenological two-state approximation omits collision trajectories,
% other relaxation mechanisms, and the remaining fluorines; even its
% fastest rates are illustrative, not an experimental liquid prediction.
% The full laboratory-frame liquid calculation disables zero-track
% elimination. Each absorptive spectrum is divided by its own maximum
% to compare line shapes, not absolute peak intensities.
%
% Calculation time: minutes.
%
% kjfritz@sandia.gov
% ilya.kuprov@weizmann.ac.il

function uf6_collisions()

% Magnetic field and representative pairs in the two chemical species
sys.magnet=0.1;
sys.isotopes={'19F','235U','19F','235U'};
inter.chem.parts={[1 2],[3 4]};

% Scalar couplings specified once per pair, in Hz
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=213;
inter.coupling.scalar{3,4}=213;

% Symmetric and collision-distorted uranium quadrupolar tensors
inter.coupling.matrix{2,2}=eeqq2nqi(0,0,7/2,[0 0 0]);
inter.coupling.matrix{4,4}=eeqq2nqi(400e6,0,7/2,[0 0 0]);

% Full laboratory-frame Redfield relaxation with fixed correlation times
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='ignore';
inter.tau_c={1e-14,1e-14};

% Complete basis and explicit numerical options
bas.formalism='sphten-liouv';
bas.approximation='none';
sys.disable={'zte'};
sys.parallel={'processes',4};

% Forward collision rates and fixed reverse rate, in Hz
collision_rates=[1e10 3e10 1e11 3e11 1e12 3e12 1e13 3e13];
reverse_rate=1e12;

% Laboratory frequency window about the unshifted fluorine resonance
parameters.spins={'19F'};
parameters.offset=0;
parameters.decouple={};
freq_centre=spin('19F')*sys.magnet/(2*pi);
parameters.sweep=freq_centre+[-1200 1200];
parameters.npoints=4801;
offset_axis=linspace(-1200,1200,parameters.npoints);
spectra=zeros(numel(collision_rates),parameters.npoints);

% Sweep the forward rate without changing the distorted-state lifetime
for n=1:numel(collision_rates)

    % Column-conserving kinetics and normalised stationary populations
    forward_rate=collision_rates(n);
    inter.chem.rates=[-forward_rate reverse_rate;...
                      forward_rate -reverse_rate];
    inter.chem.concs=[reverse_rate forward_rate]/(forward_rate+reverse_rate);

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Chemical-population-weighted excitation and unweighted detection
    parameters.rho0=state(spin_system,'L+','19F','chem');
    parameters.coil=state(spin_system,'L+','19F');

    % Laboratory-frame frequency-domain spectrum
    spectrum=liquid(spin_system,@slowpass,parameters,'labframe');
    spectra(n,:)=real(spectrum).';

end

% Unit-height normalisation for line-shape comparison
spectra=spectra./max(spectra,[],2);

% Waterfall with a logarithmic collision-rate coordinate
kfigure(); waterfall(offset_axis,log10(collision_rates),spectra,'EdgeColor','k');
kxlabel('$^{19}$F frequency offset, Hz');
kylabel('$\log_{10}(k_1/\mathrm{s}^{-1})$');
kzlabel('intensity / spectrum maximum');
kgrid; box on; axis tight;
view(-25,45);

end


