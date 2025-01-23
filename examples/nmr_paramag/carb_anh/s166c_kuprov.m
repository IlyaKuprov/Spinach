% Distributed fit for the S166C mutant dataset for human carbonic anhydrase
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

function s166c_kuprov()

% Load experimental data
load('s166c_expt.mat','expt_pcs','xyz','xyz_all');

% Load susceptibility tensor
load('s166c_chi_eff.mat','chi'); %#ok<*NODEF>

% Set inverse problem parameters
parameters.equation='kuprov';
parameters.plot={'diagnostics','density',...
                 'molecule','tightzoom','box'};
parameters.box_cent=[-14.0 -3.6 -11.0];
parameters.box_size=[ 20.0 20.0  20.0];
parameters.margins=50*ones(1,6);
parameters.confine=[3.0 12.0];
parameters.sharpen=1.0;
parameters.xyz_all=xyz_all;
parameters.expt_pcs=expt_pcs;
parameters.xyz=xyz;
parameters.chi=chi;
parameters.gpu=1;

% Solve and refine the grid
for n=[64 128 256]
    [source_cube,ranges]=ipcs(parameters,n,0.21);
    parameters.guess=source_cube;
end

% Get the new susceptibility tensor
[chi,~]=chi_eff(source_cube,ranges,xyz,expt_pcs);
disp('Effective susceptibility tensor: '); disp(chi); 

end

