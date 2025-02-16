% Figure 2 from the paper by Rasulov and Kuprov:
%
%       https://arxiv.org/abs/2502.02198
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function distortions_figure_2()

% Get E1000B pulse from Spinach
waveform=vg_pulse('E1000B',1000,0.001);
time_axis=1e3*linspace(0,0.001,1000);

% Apply amplifier compression
tanh_waveform=amp_tanh([waveform
                        zeros(size(waveform))],3e4);
tanh_waveform=tanh_waveform(1,:);
root_waveform=amp_root([waveform
                        zeros(size(waveform))],3e4,10);
root_waveform=root_waveform(1,:);

% Plot waveforms
figure(); hold on;
plot(time_axis,waveform,'LineWidth',1.5); 
plot(time_axis,tanh_waveform,'LineWidth',1.5);
plot(time_axis,root_waveform,'LineWidth',1.5);

% Plot saturation levels
plot([0 1 NaN 1 0],[3e4 3e4 NaN -3e4 -3e4],'k--');

% Annotate the plot
kxlabel('time, ms');
kylabel('nutation frequency, rad/s');
klegend({'input','output, tanh','output, $s=10$',...
         'saturation level'},'Location','SouthWest');
kgrid; xlim tight; box on; drawnow;

end

