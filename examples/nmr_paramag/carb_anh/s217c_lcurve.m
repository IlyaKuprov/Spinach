% L-curves for the S217C mutant dataset for human carbonic anhydrase
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

function s217c_lcurve()

% Load experimental data
load('s217c_expt.mat','expt_pcs','xyz','xyz_all');

% Load susceptibility tensor
load('s217c_chi_eff.mat','chi'); %#ok<*NODEF>

% Solver parameters
parameters.plot={};
parameters.equation='kuprov';
parameters.box_cent=[-21.8 -18.4  20.2];
parameters.box_size=[ 50.0  50.0  50.0];
parameters.margins=50*ones(1,6);
parameters.confine=[2.0 12.0];
parameters.sharpen=0.0;
parameters.xyz_all=xyz_all;
parameters.expt_pcs=expt_pcs;
parameters.xyz=xyz;
parameters.chi=chi;
parameters.gpu=1;

% Regularisation parameter array
lam=10.^linspace(-2.0,2.0,30);

% Result arrays
err=zeros(1,30); reg=zeros(1,30);

% Run a parallel loop
parfor n=1:30 %#ok<*PFOUS,*PFBNS>
    [~,~,~,err(n),~,reg(n)]=ipcs(parameters,64,lam(n));
    reg(n)=reg(n)/lam(n);
end

% L-curve analysis
kfigure(); S=lcurve(lam,err,reg./lam,'log'); drawnow;
disp(['Suggested smoothing parameter: ' num2str(S)]);

end

