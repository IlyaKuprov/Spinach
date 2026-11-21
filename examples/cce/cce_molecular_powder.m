% Recomputes the OX063 methylene-proton APPA curve in Fig. 9 of Jeschke
% et al., doi:10.5194/mr-6-15-2025. The complete two-proton expansion
% is evaluated for each orientation before powder averaging. No solvent
% spins are included in the simulation.
%
% ilya.kuprov@weizmann.ac.il

function cce_molecular_powder()

% OX063 methylene-proton inputs from doi:10.5281/zenodo.13850793
folder=fileparts(mfilename('fullpath'));
protons=readmatrix(fullfile(folder,'trityl_ox063_protons.csv'));
measured=readmatrix(fullfile(folder,'trityl_ox063_echo.csv'));
coords=protons(:,1:3);
hyperfine=reshape(protons(:,4:12).',3,3,[]);

% Evaluate the pair expansion before averaging the powder
times=linspace(0,40e-6,2048).';
tunnel=zeros(size(coords,1));
echo=trityl_appa(coords,hyperfine,tunnel,times);

% Compare with the measured trace in deuterated water/glycerol
kfigure();
plot(1e6*times,echo,'b-',measured(:,1),measured(:,2),'k.');
xlim([0 min(40,max(measured(:,1)))]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA simulation','Archived experimental echo');
kgrid();

end


