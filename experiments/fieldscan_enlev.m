% Plots a user-specified number of the lowest energy levels of
% the system as a function of the applied magnetic field. The 
% energies are obtained using the Arnoldi method. Syntax:
%
%         fieldscan_enlev(spin_system,parameters)
%
% Parameters:
%
%   parameters.fields  - two-element vector in Tesla,
%                        ordered as [from to]
%
%   parameters.npoints - number of points in the scan
%
%   parameters.orientation - system orientation, three-
%                            element vector containing
%                            Euler angles in radians,
%                            ordered as [alp bet gam]
%
%   parameters.nstates - number of lowest energy states
%                        to solve for
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fieldscan_enlev.m>

function fieldscan_enlev(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Generate the field grid
B0=linspace(parameters.fields(1),...
            parameters.fields(2),...
            parameters.npoints);

% Zeeman Hamiltonian at the current orientation
[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));
Hz=Hz+orientation(Qz,parameters.orientation);

% Coupling Hamiltonian at the current orientation
[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));
Hc=Hc+orientation(Qc,parameters.orientation);

% Loop over the fields
parfor n=1:numel(B0) %#ok<*PFBNS>
    
    % Build the Hamiltonian
    H=full(B0(n)*Hz+Hc);
    
    % Get level energies
    E(:,n)=eigs(H,parameters.nstates,'sr') 
    
    % Sort level energies
    E(:,n)=sort(real(E(:,n)));
    
end

% Convert the energies to cm^-1
E=hz2icm(E/(2*pi));

% Plot the energy levels
figure(); plot(B0,E,'r-'); kgrid;
kxlabel('Magnetic field, Tesla');
kylabel('Energy, cm$^{-1}$');

end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('this function is only available in zeeman-hilb formalism.');
end
if ~isfield(parameters,'fields')
    error('magnetic field ranges must be specified in parameters.fields variable.');
end
if (~isnumeric(parameters.fields))||(~isreal(parameters.fields))||...
   (~isrow(parameters.fields))||(numel(parameters.fields)~=2)||...
   (parameters.fields(1)>=parameters.fields(2))
    error('parameters.fields must be a row vector with two elements in ascending order.');
end
if spin_system.inter.magnet~=1.0
    error('sys.magnet must be set to exactly 1 Tesla for this function to work.');
end
if ~isfield(parameters,'npoints')
    error('the number of points must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (numel(parameters.npoints)~=1)||(mod(parameters.npoints,1)~=0)||...
   (parameters.npoints<1)
    error('parameters.npoints must be a positive real integer.');
end
if ~isfield(parameters,'orientation')
    error('system orientation must be specified in parameters.orientation field.');
end
if (~isnumeric(parameters.orientation))||(~isreal(parameters.orientation))||...
   (~isrow(parameters.orientation))||(numel(parameters.orientation)~=3)
    error('parameters.field must be a row vector with three elements.');
end
end

% If you take hyphens seriously, you will surely go mad.
%
% Oxford University Press style manual

