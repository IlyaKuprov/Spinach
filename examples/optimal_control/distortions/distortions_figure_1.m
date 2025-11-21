% Figure 1 from the paper by Rasulov and Kuprov:
%
%       https://arxiv.org/abs/2502.02198
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function distortions_figure_1()

% Pulse sequence and its discretisation
fapt={[0  pi/10e-6  0      0.0e-6    5.0e-6];
      [0  0         0      5.0e-6    5.5e-6];
      [0  pi/10e-6  pi/2   5.5e-6   16.0e-6];
      [0  0         0      16.0e-6  16.5e-6];
      [0  pi/10e-6  0      16.5e-6  20.0e-6]};
time_grid=linspace(-5e-6,40e-6,1000);
wave=fapt2sfo(fapt,time_grid);

% Convert waveform to mT
wave=1e3*wave/spin('1H');

% Apply a cascade of two single-pole filters
wave_spf=spf(wave,    0.9*exp(-1i*0.05));
wave_spf=spf(wave_spf,0.9*exp(-1i*0.05));

% Original vs second-order low-pass filter
kfigure(); scale_figure([1 1.6]);
subplot(3,1,1); hold on; box on;
plot(1e6*time_grid,wave(1,:),'r-','LineWidth',1);
plot(1e6*time_grid,wave(2,:),'b-','LineWidth',1);
plot(1e6*time_grid,wave_spf(1,:),'r:','LineWidth',1);
plot(1e6*time_grid,wave_spf(2,:),'b:','LineWidth',1);
kgrid; xlim tight; ylim padded;
klegend({'input, in-phase','input, quadr.',...
         'output, in-phase','output, quadr.'},...
         'Location','NorthEast');
kxlabel('time, $\mu$s'); kylabel('$B_1$, mT');

% Apply a cascade of three single-zero filters
wave_szf=szf(wave,    0.1*exp(-1i*0.05));
wave_szf=szf(wave_szf,0.1*exp(-1i*0.05));
wave_szf=szf(wave_szf,0.1*exp(-1i*0.05));

% Original vs third-order high-pass filter
subplot(3,1,2); hold on; box on;
plot(1e6*time_grid,wave(1,:),'r-','LineWidth',1);
plot(1e6*time_grid,wave(2,:),'b-','LineWidth',1);
plot(1e6*time_grid,wave_szf(1,:),'r:','LineWidth',1);
plot(1e6*time_grid,wave_szf(2,:),'b:','LineWidth',1);
kgrid; xlim tight; ylim padded;
klegend({'input, in-phase','input, quadr.',...
         'output, in-phase','output, quadr.'},...
         'Location','NorthEast');
kxlabel('time, $\mu$s'); kylabel('$B_1$, mT');

% Apply an RLC filter
wave_spf=spf(wave,    0.9*exp(-1i*0.5));
wave_spf=spf(wave_spf,0.9*exp(+1i*0.5));

% Original vs RLC filter
subplot(3,1,3); hold on; box on;
plot(1e6*time_grid,wave(1,:),'r-','LineWidth',1);
plot(1e6*time_grid,wave(2,:),'b-','LineWidth',1);
plot(1e6*time_grid,wave_spf(1,:),'r:','LineWidth',1);
plot(1e6*time_grid,wave_spf(2,:),'b:','LineWidth',1);
kgrid; xlim tight; ylim padded;
klegend({'input, in-phase','input, quadr.',...
         'output, in-phase','output, quadr.'},...
         'Location','NorthEast');
kxlabel('time, $\mu$s'); kylabel('$B_1$, mT');

end

