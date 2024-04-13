% Inverse problem for the unpaired electron density 
% distribution. Experimental data from Gottfried Ott-
% ing (Australian National University).
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk

function tm_1igv_distr_fit()

% Load the pdb file
pdb=pdbread('1igv_processed.pdb');

% Load experimental data
load('tm_1igv_pcs.mat','x','y','z','expt_pcs'); %#ok<*NODEF>

% Load susceptibility tensor
load('tm_1igv_chi_eff.mat','chi');
 
% Inverse solver parameters
parameters.equation='kuprov';
parameters.plot={'diagnostics','density',...
                 'molecule','tightzoom','box'};
parameters.box_cent=[3.5  17.0  16.1];
parameters.box_size=[7.0   7.0   7.0];
parameters.xyz_all=[[pdb.Model.Atom(:).X]'...
                    [pdb.Model.Atom(:).Y]'...
                    [pdb.Model.Atom(:).Z]'];
parameters.margins=50*ones(1,6);
parameters.confine=[1.0 3.0];
parameters.sharpen=2e3;
parameters.expt_pcs=expt_pcs;
parameters.xyz=[x y z];
parameters.chi=chi;
parameters.gpu=1;

% Iteratively refine the grid
for n=[64 128 256 384]
    [source_cube,ranges]=ipcs(parameters,n,0.34);
    parameters.guess=source_cube;
end

% Get the new susceptibility tensor
[chi,~]=chi_eff(source_cube,ranges,[x y z],expt_pcs);
disp('Effective susceptibility tensor: '); 
disp(chi); save('tm_1igv_chi_eff.mat','chi');

end

