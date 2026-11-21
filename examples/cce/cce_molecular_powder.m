% Recomputes the OX063 methylene-proton APPA curve in Fig. 9 of Jeschke
% et al., doi:10.5194/mr-6-15-2025. The complete two-proton expansion
% is evaluated for each orientation before powder averaging. No solvent
% spins are included in the simulation.
%
% Syntax: [answer,diagnostics]=cce_molecular_powder()
%
% answer is the 2048-by-1 normalised Hahn echo. diagnostics contains the
% total echo times in seconds and the powder quadrature and node traces.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_molecular_powder()

% OX063 methylene-proton inputs from doi:10.5281/zenodo.13850793
data=load(fullfile(fileparts(mfilename('fullpath')),'trityl_ox063.mat'));

% Evaluate the pair expansion before averaging the powder
times=linspace(0,40e-6,2048).';
tunnel=zeros(size(data.coords,1));
[answer,diagnostics]=trityl_appa(data.coords,data.hyperfine,tunnel,times);

% Compare with the measured trace in deuterated water/glycerol
kfigure();
plot(1e6*times,answer,'b-',data.exp_time,data.exp_echo,'k.');
xlim([0 min(40,max(data.exp_time))]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA','Experiment');
kgrid();

end


