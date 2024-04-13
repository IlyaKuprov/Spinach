% Multipolar fit for the S50C mutant dataset for human carbonic anhydrase
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
% i.kuprov@soton.ac.uk

function s50c_mult()

% Load experimental data
load('s50c_expt.mat','expt_pcs','xyz'); %#ok<*NODEF>

% Solve the inverse problem
[mxyz,chi,~,pred_pcs]=ilpcs(xyz,expt_pcs,[0 1 2],[-27.0  13.0  18.0]);

% Plot experimental vs predicted PCS
figure(); plot(expt_pcs,pred_pcs,'bo'); hold on; kgrid;
plot([min(expt_pcs) max(expt_pcs)],[min(expt_pcs) max(expt_pcs)],'r-');
xlabel('Experimental PCS, ppm'); ylabel('Predicted PCS, ppm');

% Report and save the parameters
disp('Susceptibility tensor:'); disp(chi);
disp('Magnetic multipole centre location:'); disp(mxyz);

end

