% Z magnetization of the sample as a function of magnetic field in a
% finite-speed magnetic field sweep experiment. Syntax:
%
%             magn=fieldscan_magn(spin_system,parameters)
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
%   parameters.sweep_time - sweep time, seconds
%
%   parameters.nstates - (optional) number of lowest energy 
%                        states to use for the effective 
%                        Hamiltonian in the time domain
%
% Outputs:
%
%   fields - magnetic fields in Tesla at each point in time
%
%   z_magn - total sample magnetisation at each point in time
%
% Note: this function requires Hilbert space formalism.
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fieldscan_magn.m>

function [fields,z_magn]=fieldscan_magn(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Get the rotation matrix
R=euler2dcm(parameters.orientation);

% Get total magnetic moment operator
mz=sparse(0);
for n=1:spin_system.comp.nspins
    
    % Get the rotated g-tensor
    g=R*gtensorof(spin_system,n)*R';
   
    % Biild elementary operators
    Sp=operator(spin_system,{'L+'},{n}); Sm=Sp';
    Sx=(Sp+Sm)/2; Sy=(Sp-Sm)/2i; Sz=-1i*(Sx*Sy-Sy*Sx);
    
    % Build magnetic moment operator
    mz=mz-g(3,1)*Sx-g(3,2)*Sy-g(3,3)*Sz;
    
end

% Generate the field grid
fields=linspace(parameters.fields(1),...
                parameters.fields(2),...
                parameters.npoints);

% Zeeman Hamiltonian
[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));

% Coupling Hamiltonian
[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));

% Thermal equilibrium at the starting field
rho=equilibrium(spin_system,parameters.fields(1)*Hz+Hc,...
                            parameters.fields(1)*Qz+Qc,...
                            parameters.orientation);
rho=rho/trace(rho);

% Preallocate the answer
z_magn=zeros(size(fields),'like',1i);

% Hamiltonians at the current orientation
Hz=Hz+orientation(Qz,parameters.orientation);
Hc=Hc+orientation(Qc,parameters.orientation);

% Project into the active space
if isfield(parameters,'nstates')
    [V,~]=eigs(mean(fields)*Hz+Hc,parameters.nstates,'sr');
    Hz=V'*Hz*V; Hc=V'*Hc*V; rho=V'*rho*V; mz=V'*mz*V;
end

% Time stepping
dt=parameters.sweep_time/(parameters.npoints-1);

% Evolution and acquisition
for n=1:parameters.npoints
    
    % Compute the observable
    z_magn(n)=hdot(rho,mz);
    
    % Compute the propagator
    P=propagator(spin_system,fields(n)*Hz+Hc,dt);
    
    % Apply the propagator
    rho=P*rho*P';

end

% Discard the imaginaries
z_magn=real(z_magn);

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
if ~isfield(parameters,'npoints')
    error('the number of points must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (numel(parameters.npoints)~=1)||(mod(parameters.npoints,1)~=0)||...
   (parameters.npoints<2)
    error('parameters.npoints must be a real integer greater than 1.');
end
if ~isfield(parameters,'orientation')
    error('system orientation must be specified in parameters.orientation field.');
end
if (~isnumeric(parameters.orientation))||(~isreal(parameters.orientation))||...
   (~isrow(parameters.orientation))||(numel(parameters.orientation)~=3)
    error('parameters.field must be a row vector with three elements.');
end
if ~isfield(parameters,'sweep_time')
    error('sweep time must be specified in parameters.sweep_time field.');
end
if (~isnumeric(parameters.sweep_time))||(~isreal(parameters.sweep_time))||...
   (numel(parameters.sweep_time)~=1)||(parameters.sweep_time<=0)
    error('parameters.sweep_time must be a positive real number.');
end
end

% Those who are determined to be 'offended' will discover a provocation
% somewhere. We cannot possibly adjust enough to please the fanatics, and
% it is degrading to make the attempt.
%
% Christopher Hitchens

