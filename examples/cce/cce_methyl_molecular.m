% Recomputes the 36-proton Finland trityl core APPA curve in Fig. 7 of
% Jeschke et al., doi:10.5194/mr-6-15-2025. Methyl tunnelling enters as
% coherent intramethyl exchange, not as stochastic jumps. The measured
% echo is overlaid with the scale factor used by the authors.
%
% ilya.kuprov@weizmann.ac.il

function cce_methyl_molecular()

% Finland trityl core inputs from doi:10.5281/zenodo.13850793
folder=fileparts(mfilename('fullpath'));
protons=readmatrix(fullfile(folder,'trityl_finland_protons.csv'));
measured=readmatrix(fullfile(folder,'trityl_finland_echo.csv'));
coords=protons(:,1:3);
hyperfine=reshape(protons(:,4:12).',3,3,[]);

% Three protons per methyl group share a tunnel splitting (kHz)
splittings=[54.5 4.8 4.8 54.5 54.5 4.8 54.5 4.8 4.8 54.5 54.5 4.8];
tunnel=kron(diag(2*pi*1e3*splittings),ones(3)-eye(3));

% Evaluate the pair expansion before averaging the powder
times=linspace(0,40e-6,2048).';
echo=trityl_appa(coords,hyperfine,tunnel,times);

% Compare with the measured trace in deuterated o-terphenyl
kfigure();
plot(1e6*times,echo,'b-',measured(:,1),measured(:,2),'k.');
xlim([0 max(measured(:,1))]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA simulation (36 protons)','Archived experimental echo');
kgrid();

end


