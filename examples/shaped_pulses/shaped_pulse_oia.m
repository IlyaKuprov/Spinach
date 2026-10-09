% Offset-independent adiabaticity inversion pulses from Table 1 of
% Tannus and Garwood (JMR A 120, 133 (1996)): amplitude and frequ-
% ency modulation functions, and the inversion profiles of a single
% proton, for a 50 kHz sweep in 2 ms. Reproduces Figure 1 of the
% paper on a +/-30 kHz offset grid.
%
% Calculation time: 16 seconds on 11 workers, most of
%                   it the parallel pool startup
%
% ilya.kuprov@weizmann.ac.il

function shaped_pulse_oia()

% Single proton
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Operators and the initial state
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');
rho0=state(spin_system,'Lz','1H');
rho0=rho0/norm(rho0,2);

% Amplitude functions from Table 1 of the paper
am_funs={@(tau)1./(1+99*tau.^2),...
         @(tau)sech(asech(0.01)*tau),...
         @(tau)exp(-log(100)*tau.^2),...
         @(tau)(1+cos(pi*tau))/2,...
         @(tau)sech(asech(0.01)*tau.^8),...
         @(tau)1-abs(sin(pi*tau/2)).^40};
am_names={'Lorentz','HS','Gauss','Hanning','HS8','Sin40'};

% Pulse parameters from the paper
npts=1000; dur=2e-3; bwidth=50e3;

% Offset grid, Hz
offsets=linspace(-1.2,1.2,121)*bwidth/2;

% Loop over the amplitude functions
kfigure(); scale_figure([2.00 0.75]);
for n=1:numel(am_funs)

    % Make the pulse
    [Cx,Cy,durs,~,amps,~,frqs]=oia_pulse(npts,dur,bwidth,am_funs{n});
    time_axis=1e3*cumsum(durs);

    % Inversion profile
    mz=zeros(size(offsets));
    parfor k=1:numel(offsets)
        rho=shaped_pulse_xy(spin_system,2*pi*offsets(k)*Lz,{Lx,Ly},{Cx,Cy},durs,rho0,'expv-pwc');
        mz(k)=real(rho0'*rho);
    end

    % Plotting
    subplot(1,3,1); plot(time_axis,amps/(2*pi*1e3)); hold on;
    subplot(1,3,2); plot(time_axis,frqs/1e3); hold on;
    subplot(1,3,3); plot(offsets/1e3,mz); hold on;

end

% Axis labels and legends
subplot(1,3,1); kxlabel('time, ms'); kylabel('amplitude, kHz'); ylim padded;
ktitle('amplitude modulation'); klegend(am_names,'Location','northeast'); kgrid; xlim tight;
subplot(1,3,2); kxlabel('time, ms'); kylabel('frequency, kHz'); ylim padded;
ktitle('frequency modulation'); klegend(am_names,'Location','northwest'); kgrid; xlim tight;
subplot(1,3,3); kxlabel('offset, kHz'); kylabel('Mz/M0'); ylim padded;
ktitle('inversion profile'); klegend(am_names,'Location','north'); kgrid; xlim tight;

end


