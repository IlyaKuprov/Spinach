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
% stepped through by the pulse sequence, which averages over the rotor
% phase at the start of the sequence and keeps the electron coherence
% pathway (-1 after the first pulse, +1 after the second) in place of
% the phase cycle. All P1 centre tensors are axial and coaxial, a two-
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

% Integrated two-pulse echo as a function of the carrier offset
function echo=echo_sweep(spin_system,parameters,H,~,~)

% Carrier offsets across the sweep
offsets=ft_axis(0,parameters.sweep,parameters.npoints);

% Pulse and offset operators
sx=operator(spin_system,'Lx','E');
sz=operator(spin_system,'Lz','E');

% Step counts of the pulses, the delay, and the echo window
pulse_steps=round(parameters.pulse_dur/parameters.timestep);
delay_steps=round(parameters.tau/parameters.timestep);
echo_steps=round(parameters.echo_win/parameters.timestep);
nsteps=2*pulse_steps+delay_steps+echo_steps;

% Rotor stack advance at the middle of each time step
stack_shift=round(parameters.rate*parameters.timestep*((1:nsteps)-1/2)*parameters.spc_dim);

% Rotor stack indices at the start of the sequence
start_idx=floor((0:(parameters.nphases-1))*parameters.spc_dim/parameters.nphases);

% Free evolution propagators at every rotor phase
p_free=cell(parameters.spc_dim,1);
for n=1:parameters.spc_dim
    p_free{n}=propagator(spin_system,H{n},parameters.timestep);
end

% Preallocate the answer
echo=zeros(parameters.npoints,1);

% Loop over carrier offsets
for k=1:parameters.npoints

    % Carrier offset propagator, the rotor stack commutes with sz
    p_off=propagator(spin_system,2*pi*offsets(k)*sz,parameters.timestep);

    % Pulse propagators at every rotor phase
    p_pulse=cell(parameters.spc_dim,1);
    for n=1:parameters.spc_dim
        p_pulse{n}=propagator(spin_system,H{n}+2*pi*offsets(k)*sz+...
                              2*pi*parameters.pulse_frq*sx,parameters.timestep);
    end

    % Loop over rotor phases at the start of the sequence
    for j=1:parameters.nphases

        % Rotor stack indices at each time step
        idx=mod(start_idx(j)+stack_shift,parameters.spc_dim)+1;

        % First pulse
        rho=parameters.rho0;
        for s=1:pulse_steps
            rho=p_pulse{idx(s)}*rho*p_pulse{idx(s)}';
        end

        % Select the -1 coherence on the electron
        rho=coherence(spin_system,rho,{{'E',-1}});

        % Interpulse delay
        for s=(pulse_steps+1):(pulse_steps+delay_steps)
            rho=p_off*p_free{idx(s)}*rho*p_free{idx(s)}'*p_off';
        end

        % Second pulse
        for s=(pulse_steps+delay_steps+1):(2*pulse_steps+delay_steps)
            rho=p_pulse{idx(s)}*rho*p_pulse{idx(s)}';
        end

        % Select the +1 coherence on the electron
        rho=coherence(spin_system,rho,{{'E',+1}});

        % Integrate the signal over the echo window
        for s=(2*pulse_steps+delay_steps+1):nsteps
            rho=p_off*p_free{idx(s)}*rho*p_free{idx(s)}'*p_off';
            echo(k)=echo(k)+trace(parameters.coil'*rho);
        end

    end

end

end

