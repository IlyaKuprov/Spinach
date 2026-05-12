% Builds a small quiet Spinach spin system for tests. Syntax:
%
%                    spin_system=test_spin_system(sys,inter,bas)
%
% Parameters:
%
%     sys          - Spinach system specification
%
%     inter        - Spinach interaction specification
%
%     bas          - Spinach basis specification
%
% Outputs:
%
%     spin_system  - Spinach spin system object
%
% ilya.kuprov@weizmann.ac.il

function spin_system=test_spin_system(sys,inter,bas)

% Apply quiet settings used by regression tests
sys.output='hush';
if isfield(sys,'disable')
    sys.disable=unique([sys.disable {'hygiene'}]);
else
    sys.disable={'hygiene'};
end
sys.parallel={'local',1};
sys.parprops={};

% Build the Spinach object and basis
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

end
