% Recomputes the 36-proton Finland trityl core APPA curve in Fig. 7 of
% Jeschke et al., doi:10.5194/mr-6-15-2025. Methyl tunnelling enters as
% coherent intramethyl exchange, not as stochastic jumps. The measured
% echo is overlaid with the scale factor used by the authors.
%
% ilya.kuprov@weizmann.ac.il

function cce_methyl_molecular()

% Finland trityl core inputs from doi:10.5281/zenodo.13850793
data=load(fullfile(fileparts(mfilename('fullpath')),'trityl_finland.mat'));

% Evaluate the pair expansion before averaging the powder
times=linspace(0,40e-6,2048).';
echo=trityl_appa(data.coords,data.hyperfine,data.tunnel,times);

% Compare with the measured trace in deuterated o-terphenyl
kfigure();
plot(1e6*times,echo,'b-',data.exp_time,data.exp_echo,'k.');
xlim([0 max(data.exp_time)]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA simulation (36 protons)','Archived experimental echo');
kgrid();

end


