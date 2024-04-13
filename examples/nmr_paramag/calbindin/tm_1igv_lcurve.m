% Inverse problem for the unpaired electron density 
% distribution. Experimental data from Gottfried Ott-
% ing (Australian National University).
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk

function tm_1igv_lcurve()

% Load the pdb file
pdb=pdbread('1igv_processed.pdb');

% Load experimental data
load('tm_1igv_pcs.mat','x','y','z','expt_pcs'); %#ok<*NODEF>

% Load susceptibility tensor
load('tm_1igv_chi_eff.mat','chi');

% Inverse solver parameters
parameters.plot={};
parameters.equation='kuprov';
parameters.box_cent=[3.5  17.0  16.1];
parameters.box_size=[7.0   7.0   7.0];
parameters.xyz_all=[[pdb.Model.Atom(:).X]'...
                    [pdb.Model.Atom(:).Y]'...
                    [pdb.Model.Atom(:).Z]'];
parameters.margins=50*ones(1,6);
parameters.confine=[1.0 3.0];
parameters.sharpen=0.0;
parameters.expt_pcs=expt_pcs;
parameters.xyz=[x y z];
parameters.chi=chi;
parameters.gpu=1;
 
% Regularisation parameter array
lam=10.^linspace(-1.5,1.0,15);

% Result arrays
err=zeros(1,15); reg=zeros(1,15);

% Run a parallel loop
parfor n=1:15 %#ok<*PFOUS,*PFBNS>
    [~,~,~,err(n),~,reg(n)]=ipcs(parameters,128,lam(n));
    reg(n)=reg(n)/lam(n);
end

% L-curve analysis
S=lcurve(lam,err,reg,'log'); drawnow;
disp(['Suggested smoothing parameter: ' num2str(S)]); 
          
end

