% Point fit for the S220C mutant dataset for human carbonic anhydrase
% II. The system and the method are described in:
%
%                   http://dx.doi.org/10.1039/c6sc03736d
%
% A step-by-step tutorial is available here:
%
% http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis
%
% e.suturina@soton.ac.uk
% daniel.haeussinger@unibas.ch
% kaspar.zimmermann@unibas.ch
% maxim.yulikov@phys.chem.ethz.ch
% luca.garbuio@psi.ch
% gunnar.jeschke@phys.chem.ethz.ch
% ilya.kuprov@weizmann.ac.il

function s220c_point()

% Load experimental data
load('s220c_expt.mat','expt_pcs','xyz');

% Solve the inverse problem
[mxyz,chi,pred_pcs]=ippcs(xyz,[-14 -26 4],expt_pcs);

% Plot experimental vs predicted PCS
figure(); plot(expt_pcs,pred_pcs,'bo'); hold on; kgrid;
plot([min(expt_pcs) max(expt_pcs)],[min(expt_pcs) max(expt_pcs)],'r-');
kxlabel('Experimental PCS, ppm'); kylabel('Predicted PCS, ppm');
xlim([min([expt_pcs; pred_pcs]) max([expt_pcs; pred_pcs])]);
ylim([min([expt_pcs; pred_pcs]) max([expt_pcs; pred_pcs])]);

% Report and save the parameters
disp('Susceptibility tensor:'); disp(chi);
disp('Point electron location:'); disp(mxyz);

end

