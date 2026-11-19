% Hahn echo with a frequency-selected bath pump and finite bandwidth.
% The pump selects physical bath site 2 in a 200 kHz spectral window.
% The deliberately small model is not a published material reproduction.
% Syntax:
%
%             [answer,diagnostics]=cce_deer_bandwidth()
%
% Returns the normalised complex signal column in answer and the calculation
% diagnostics structure in diagnostics, with one entry per time point.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_deer_bandwidth()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','E','E'};
sys.output='hush';
sys.parallel={'local',1};

% Toy coordinates in Angstrom and zero electron g tensors
inter.coordinates={[0 0 0];[25 0 0];[0 35 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit electron-electron coupling tensors in Hz
inter.coupling.matrix{1,2}=diag([0 0 20000]);
inter.coupling.matrix{1,3}=diag([0 0 12000]);
inter.coupling.matrix{2,3}=diag([200 200 -400]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Experiment parameters and retained bath order
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.times=linspace(0,200e-6,41);
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);
parameters.pump=struct('sites',2,'center',0,'halfwidth',1e5,...
                       'units','Hz','fraction',0.25,'allow_empty',true);

% Run the cluster expansion for the named bath-pumped Hahn experiment
[answer,diagnostics]=clexpn(spin_system,@bath_pumped_hahn,parameters,'labframe');

% Plot the projected complex signal
kfigure();
plot(1e6*parameters.times,real(answer),'b-',...
     1e6*parameters.times,imag(answer),'r--');
xlabel('total echo time / microseconds');
ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary');
kgrid();

end


