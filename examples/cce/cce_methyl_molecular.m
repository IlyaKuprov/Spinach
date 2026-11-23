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

% High-field dipolar geometry and electron-proton tensors in Spinach input units
sys.magnet=1.2414;
sys.isotopes=[{'E'} repmat({'1H'},1,size(protons,1))];
sys.output='hush';
sys.tols.prox_cutoff=0;
sys.tols.cce_atol=0;
sys.tols.cce_mem=2^27;
inter.coordinates=[{[]}; num2cell(protons(:,1:3),2)];
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
for k=1:size(protons,1)
    inter.coupling.matrix{1,k+1}=reshape(protons(k,4:12),3,3)/(2*pi);
end

% Three protons per methyl group share a tunnel splitting (kHz)
splittings=[54.5 4.8 4.8 54.5 54.5 4.8 54.5 4.8 4.8 54.5 54.5 4.8];
inter.coupling.scalar=cell(numel(sys.isotopes),numel(sys.isotopes));
for k=1:numel(splittings)
    sites=3*k+(-1:1);
    inter.coupling.scalar(sites(1),sites(2:3))={-2e3*splittings(k)/3};
    inter.coupling.scalar{sites(2),sites(3)}=-2e3*splittings(k)/3;
end

% Spinach housekeeping
spin_system=create(sys,inter);

% Complete secular pair expansion at every powder node
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[1/2 -1/2]);
parameters.times=linspace(0,40e-6,2048);
parameters.grid='icos_2ang_642pts';
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.backend='analytical_hahn';

% Simulation
echo=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Compare with the measured trace in deuterated o-terphenyl
kfigure();
plot(1e6*parameters.times,echo,'b-',measured(:,1),measured(:,2),'k.');
xlim([0 max(measured(:,1))]);
ylim([0 1]);
kxlabel('Total echo time / $\mu$s');
kylabel('Normalised Hahn echo');
klegend('APPA simulation (36 protons)','Archived experimental echo');
kgrid();

end


