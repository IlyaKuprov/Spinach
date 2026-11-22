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

% High-field dipolar geometry and electron-proton tensors in Spinach input units
sys.magnet=1.2414;
sys.isotopes=[{'E'} repmat({'1H'},1,size(protons,1))];
sys.output='hush';
sys.tols.prox_cutoff=0;
inter.coordinates=[{[]}; num2cell(protons(:,1:3),2)];
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
for k=1:size(protons,1)
    inter.coupling.matrix{1,k+1}=reshape(protons(k,4:12),3,3)/(2*pi);
end

% Evaluate the complete secular pair expansion at every powder node before averaging
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[1/2 -1/2]);
parameters.times=linspace(0,40e-6,2048);
parameters.grid='icos_2ang_642pts';
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.backend='analytical_hahn';
parameters.cce.numerics=struct('atol',0,'max_memory_bytes',2^27);
spin_system=create(sys,inter);
echo=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Compare with the measured trace in deuterated water/glycerol
kfigure();
plot(1e6*parameters.times,echo,'b-',measured(:,1),measured(:,2),'k.');
xlim([0 min(40,max(measured(:,1)))]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA simulation','Archived experimental echo');
kgrid();

end


