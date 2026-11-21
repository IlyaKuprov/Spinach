% Finite observed spin-lock drive including the coupled nuclear drift.
% Piecewise constant drive in rad/s; dense reference is in test_cce_examples.
% This small illustrative model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=gcce_drive_small()
%
% Returns one row of column-vectorised observed density per time point
% in answer, and the cluster-expansion diagnostics structure in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_drive_small()

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0;
sys.isotopes={'E','13C'};
sys.output='hush';
sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,2);
inter.coupling.matrix=cell(2);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([1100 900 2200]);

% Compile metadata, without a full-bath basis
spin_system=create(sys,inter);

% Basis formalism for retained clusters
bas.formalism='zeeman-hilb';

% Experiment parameters
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.preparation='transition';
parameters.output='density';
parameters.times=[0 linspace(20e-6,200e-6,30)];
parameters.drive=struct('spins',1,'axis','x','amplitude',4e4,...
                        'detuning',0,'duration',20e-6);
parameters.cce.bath_spins=2;
parameters.cce.approximation='generalized';
parameters.cce.order=1;
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Run the named finite-drive experiment
[answer,diagnostics]=clexpn(spin_system,@finite_drive,parameters,'labframe');

end


