% Projected Hahn echo of the original single-group P1 numerical fixture.
% Syntax: [answer,diagnostics]=p1_hahn_fixture()
% Returns the complex column on a 41-point 0 to 200 microsecond grid
% and the native clexpn audit for the complete electron/14N group.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=p1_hahn_fixture()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','E','14N'};
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-3;
sys.tols.cce_rtol=1e-3;
sys.tols.cce_mem=2^28;

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[35 0 0];[35 0 1]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=diag([2000 2000 -4000]);
inter.coupling.matrix{2,3}=diag([81 81 114])*1e6;
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;400 0 200].';

% Declare observed-Sz projection for this projected electron-bath model
for n=2:numel(sys.isotopes)

    % Remove transverse observed-spin coupling components
    if ~isempty(inter.coupling.matrix{1,n})
        inter.coupling.matrix{1,n}(:,1:2)=0;
    end

end

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Spinach housekeeping
spin_system=create(sys,inter);

% Projected transition and total echo times in seconds
parameters.observed=1;
parameters.cce.bath_spins=2:3;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,200e-6,41);
parameters.cce.order=1;
parameters.cce.approximation='projected';

% Keep the electron and its nitrogen in the same physical group
parameters.cce.expansion='partition';
parameters.cce.groups={[2 3]};
parameters.cce.partition=struct('mode','user','max_dimension',6);

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end


