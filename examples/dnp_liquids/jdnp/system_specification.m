% Parameters of the 2e1n system used for the simulations reported
% in https://doi.org/10.1039/d1cp04186j
%
% maria-grazia.concilio@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function [sys,inter,bas,parameters]=system_specification()

% Spin system
sys.isotopes={'1H','E','E'};            
               
% Nuclear chemical shift tensor
inter.zeeman.matrix{1}=diag([5 10 20]); 

% Electron g-tensors - axial along Z
inter.zeeman.matrix{2}=diag([2.0032 2.0032 2.0026]);
inter.zeeman.matrix{3}=diag([2.0032 2.0032 2.0026]);
                    
% Set coordinates           
inter.coordinates{1,1}=[-3.00 0.50  1.30];
inter.coordinates{2,1}=[ 0.00 0.00 -9.37];
inter.coordinates{3,1}=[ 0.00 0.00 +9.37];

% Empty scalar coupling array
inter.coupling.scalar=cell(3,3); 
         
% Relaxation theory
inter.relaxation={'redfield','SRFK'};
inter.equilibrium='dibari';
inter.rlx_keep='labframe';
inter.temperature=298;
inter.tau_c={500e-12};
inter.srfk_tau_c={[1.0 1e-12]};
inter.srfk_mdepth=cell(3);
inter.srfk_mdepth{2,3}=3e9;    

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Tolerance settings
sys.tols.rlx_integration=1e-10;
sys.disable={'hygiene'};
sys.output='hush';

% Reference g-factors
parameters.g_ref=2.00231930436256;
parameters.g_trityl=mean(diag(inter.zeeman.matrix{2}));

end