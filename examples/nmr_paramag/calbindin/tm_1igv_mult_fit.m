% Electron location and susceptibility tensor recovery from
% experimental PCS data using point electron model. Experi-
% mental data kindly provided by Gottfried Otting (Australi-
% an National University).
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk

function tm_1igv_mult_fit()

% Load experimental data
load('tm_1igv_pcs.mat','expt_pcs','x','y','z');

% Solve the inverse problem
[mxyz,chi,~,pred_pcs]=ilpcs([x y z],expt_pcs,[0 1 2],[-5 5 -15]);

% Plot experimental vs predicted PCS
figure(); plot(expt_pcs,pred_pcs,'bo'); hold on; kgrid;
plot([min(expt_pcs) max(expt_pcs)],[min(expt_pcs) max(expt_pcs)],'r-');
kxlabel('Experimental PCS, ppm'); kylabel('Predicted PCS, ppm');
xlim([min([expt_pcs; pred_pcs]) max([expt_pcs; pred_pcs])]);
ylim([min([expt_pcs; pred_pcs]) max([expt_pcs; pred_pcs])]);

% Report the parameters
disp('Susceptibility tensor:'); disp(chi);
disp('Point electron location:'); disp(mxyz);

end

