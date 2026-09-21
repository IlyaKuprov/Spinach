% Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework,
% a J=8 giant spin with a crystal field to twelfth spherical rank, under
% four magnetic field profiles: linear sweep, piecewise linear sweep,
% monotone cubic spline through a measured 65 T short pulse, and a sinu-
% soidal field at the clock transition frequency. Spin-phonon relaxation
% is the generalised Lindblad dissipator of Saito and Miyashita with a
% super-Ohmic phonon bath. Reproduces Figure 2 of
%
%         https://arxiv.org/abs/2609.16352
%
% with the crystal field parameters, g-factor, temperatures, spectral
% density, sweep profiles, and stair widths of that paper. As in the
% paper, the first three profiles are propagated for 1 ms (the first
% millisecond of the measured 10 ms pulse), the sinusoid for 140 ps.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function ho_pzdo4_profiles()

% Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
[ks,qs,bkq]=ho_pzdo4_params();

% Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank
coeff=cell(1,12); euler=cell(1,12);
for k=1:12
    stev=zeros(2*k+1,1); sel=(ks==k); stev(qs(sel)+k+1)=bkq(sel);
    coeff{k}=stev2sph(k,icm2hz(stev)); euler{k}=[0 0 0];
end

% Magnet must be 1 Tesla, the field is set by the sweep
sys.magnet=1.0;

% Parallel pool size
sys.parallel={'processes',4};

% J=8 giant spin, effective g-factor 1.24
sys.isotopes={'E17'};
inter.zeeman.scalar={1.24};
inter.giant.coeff={coeff};
inter.giant.euler={euler};

% Formalism and basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spin-phonon coupling operator: unit elements between adjacent m_J states
Jz=full(stevens(17,1,0)); mj=diag(Jz);
parameters.phonon_x=double(abs(mj-mj.')==1);

% Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units
parameters.phonon_alpha=2;
parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2;

% Observable: magnetic moment along Z in Bohr magnetons
parameters.coil=-1.24*Jz;

% Single crystal, crystal field frame aligned with the laboratory frame
parameters.spins={'E17'}; parameters.orientation=[0 0 0];
parameters.needs={'zeeman_op'};

% Measured 65 T pulse of the paper, ms and Tesla, 76 points, for the monotone spline profile
pulse_t=[0 0.2 0.3 0.4 0.5 0.7 0.8 1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 ...
         2.9 3 3.2 3.3 3.4 3.5 3.7 3.8 4 4.1 4.2 4.3 4.4 4.7 4.8 5 5.1 5.2 5.5 5.6 5.8 5.9 6 ...
         6.1 6.3 6.5 6.6 6.7 6.8 6.9 7 7.4 7.6 7.7 7.8 7.9 8 8.3 8.4 8.5 8.6 8.8 8.9 9 9.2 9.3 ...
         9.4 9.5 9.6 9.7 10]*1e-3;
pulse_b=[0.000000 2.136074 3.642848 5.179167 6.810767 8.966783 10.895307 12.742582 14.767864 16.347022 ...
         17.706074 19.153758 20.263158 21.145803 22.123729 23.408180 24.796775 25.978559 ...
         27.197273 28.194403 29.228464 30.308319 31.066138 31.887477 33.187439 34.191955 ...
         35.639640 37.264593 38.557907 39.716794 40.710969 42.050817 43.084877 44.478643 ...
         45.444013 46.088824 47.040160 48.254442 49.908201 51.357363 52.539147 53.189128 ...
         53.839109 55.464061 56.054953 56.930950 57.548432 58.093530 58.979867 59.984384 ...
         60.722998 61.166167 61.520702 61.918077 62.277044 62.791120 63.707002 64.120626 ...
         64.297894 64.475161 64.681974 64.829697 65.154687 65.243321 65.272865 65.302410 ...
         65.125142 65.006964 65.006964 64.829697 64.741063 64.622884 64.475161 64.297894 ...
         64.120626 63.588824];

% Four field profiles, Tesla as a function of time in seconds
profiles={@(t) 1e4*t, ...
          @(t) interp1([0 1e-6 1e-5 1e-4 1e-3],[0 0.1 1 5 10],t,'linear'), ...
          @(t) pchip(pulse_t,pulse_b,t), ...
          @(t) 0.1*sin(0.134124264765e12*t)};
temps=[2.0 2.0 2.0 0.001]; steps=[1e-8 1e-8 1e-8 1e-15]; nsteps=[1e5 1e5 1e5 140538]; nout=[100 100 100 59];
labels={'linear, 10 T/ms','piecewise linear','spline of a measured pulse','sinusoidal, 0.1 T at the clock gap'};
tscale=[1e6 1e6 1e6 1e12]; tunits={'$\mu$s','$\mu$s','$\mu$s','ps'};

% Loop over the profiles
kfigure(); scale_figure([2.0 1.6]); answers=cell(1,4);
for n=1:4

    % Spinach housekeeping at the temperature of the panel
    inter.temperature=temps(n);
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Sweep parameters of the panel
    parameters.field_prof=profiles{n};
    parameters.timestep=steps(n);
    parameters.nsteps=nsteps(n);
    parameters.nout=nout(n);

    % Run the simulation
    answers{n}=crystal(spin_system,@pulsed_field,parameters,'labframe');

    % Plot the field and the magnetisation against time
    subplot(2,2,n); yyaxis left; plot(answers{n}.t*tscale(n),answers{n}.field); kylabel('Field, Tesla');
    yyaxis right; plot(answers{n}.t*tscale(n),answers{n}.obs); kylabel('Magnetisation, $\mu_B$');
    kxlabel(['Time, ' tunits{n}]); ktitle(labels{n}); kgrid; xlim tight; drawnow;

end

% Save the curves
save('ho_pzdo4_profiles.mat','answers','labels');

end

