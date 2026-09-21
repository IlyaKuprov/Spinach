% Two-pulse echo-detected frequency-swept EPR spectra of the P1 sub-
% stitutional nitrogen defect in diamond, static and under magic ang-
% le spinning, after Figure 1a of Khamrui et al., J. Phys. Chem. Lett.
% 2026, <https://doi.org/10.1021/acs.jpclett.6c02108>: 400 ns long
% pulses with 416 kHz nutation frequency, 300 ns interpulse delay, at
% 6.9 T, static and at 10, 25, and 37 kHz MAS. The carrier is stepped
% across the spectrum and the integrated echo is recorded at each fre-
% quency, with the spectra normalised to the static one.
%
% The 14N hyperfine coupling (dipolar part 10.9 MHz) makes the two outer
% lines dephase under spinning as their resonance frequencies move during
% the sequence, whereas the central line survives. Compared to the paper's
% own simulation, the static outer edges here are weaker (0.2 of the cent-
% ral peak against 0.4) and the central line dephases less (0.8 at 37 kHz).
%
% The Hamiltonian rotor stack built by singlerot.m in Hilbert space is
% stepped through by echo_sweep.m, which averages over the rotor phase
% at the start of the sequence and keeps the electron coherence path-
% way (-1 after the first pulse, +1 after the second) in place of the
% phase cycle. All P1 centre tensors are axial and coaxial, a two-
% angle powder grid is therefore sufficient. Relaxation is not inclu-
% ded because it scales the four spectra by the same factor, which the
% normalisation removes. Slow spinning needs a high rotor rank because
% the stack must resolve the rotor phase to within one time step.
%
% Calculation time: hours on a 256-core node.
%
% ilya.kuprov@weizmann.ac.il

function mas_diamond_p1()

% P1 centre parameters
p1_params.orientation='111';
p1_params.nitrogen='14N';

% Build the spin system
[sys,inter]=diamond_p1(p1_params);

% Magnet field, central line at 193.797 GHz
sys.magnet=6.9156;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotor parameters
parameters.axis=[1 1 1];
parameters.max_rank=2700;

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.pulse_dur=400e-9;
parameters.pulse_frq=416e3;
parameters.tau=300e-9;
parameters.echo_win=1.0e-6;
parameters.timestep=5e-9;
parameters.nphases=100;
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=601;
parameters.zerofill=601;
parameters.grid='rep_2ang_400pts_sph';
parameters.axis_units='GHz-labframe';
parameters.verbose=0;

% Spinning rates
rates=[0 10e3 25e3 37e3];

% Simulation
spectra=zeros(parameters.npoints,numel(rates));
for n=1:numel(rates)
    parameters.rate=rates(n);
    spectra(:,n)=abs(singlerot(spin_system,@echo_sweep,parameters,'esr'));
end

% Normalisation to the static spectrum
spectra=spectra/max(spectra(:,1));

% Plotting
kfigure(); hold on;
for n=1:numel(rates)
    plot_1d(spin_system,spectra(:,n),parameters);
end
klegend({'static','10 kHz MAS','25 kHz MAS','37 kHz MAS'},'Location','NorthEast');
kylabel('echo intensity, a.u.');

end


