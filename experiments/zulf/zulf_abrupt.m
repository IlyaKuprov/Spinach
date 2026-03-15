% Zero-field magnetometry experiment that propagates the initial condition
% through an exponential drop in the external magnetic field and then runs
% the detection at zero field. Syntax:
%
%                fid=zulf_abrupt(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    .drop_field   - the magnetic field that the sample should 
%                    be dropped to, starting from the field spe-
%                    cified in sys.magnet, Tesla
%
%    .drop_time    - drop time, seconds
%
%    .drop_npoints - number of discretisation points in the drop
% 
%    .drop_rate    - field drop rate, Hz
%
%    .sweep        - sweep width during acquisition
%
%    .npoints      - number of points during acquisition
%
%    .detection    - 'uniaxial' to emulate common ZULF
%                    hardware, 'quadrature' for proper
%                    frequency sign discrimination
%
%    .flip_angle   - pulse flip angle in radians for
%                    protons; for other nuclei, this
%                    will be scaled by the gamma ratio
%
% Outputs:
%
%     fid          - the free induction decay detected on the
%                    internally generated gamma-weighted state
%
% Note: this function ignores the offset parameter and makes its
%       own Hamiltonian.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=zulf_abrupt.m>

function fid=zulf_abrupt(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Build the coupling Hamiltonian
H_z=hamiltonian(assume(spin_system,'labframe','zeeman'));
H_c=hamiltonian(assume(spin_system,'labframe','couplings'));

% Normalise Zeeman to 1 Tesla
H_z=H_z/spin_system.inter.magnet;

% Generate an exponential drop
drop=expdrop(spin_system.inter.magnet,parameters.drop_field,parameters.drop_time,...
             parameters.drop_npoints,parameters.drop_rate);
time_step=parameters.drop_time/parameters.drop_npoints;

% Isotropic thermal equilibrium
rho=equilibrium(spin_system);
         
% Propagate the system through the drop
for n=1:numel(drop)
    
    % Take a propagation step
    report(spin_system,['magnetic field drop step ' num2str(n) '/' num2str(numel(drop)) '...']);
    rho=evolution(spin_system,drop(n)*H_z+H_c+1i*R+1i*K,[],rho,time_step,1,'final');
    
end

% Magnetogyric ratio weights relative to 1H
weights=spin_system.inter.gammas/spin('1H');

% Get gamma-weighted detection state
coil=sparse(0);
for n=1:spin_system.comp.nspins
    coil=coil+weights(n)*state(spin_system,{'L+'},{n});
end

% Get gamma-weighted pulse operator
Sy=sparse(0);
for n=1:spin_system.comp.nspins
    Sy=Sy+weights(n)*(operator(spin_system,{'L+'},{n})-...
                      operator(spin_system,{'L-'},{n}))/2i;
end

% Apply the pulse
rho=step(spin_system,Sy,rho,parameters.flip_angle);

% Compute the digitization parameter
timestep=1/parameters.sweep;

% Run the detection at zero field
fid=evolution(spin_system,H_c+1i*R+1i*K,coil,rho,timestep,parameters.npoints-1,'observable');

% Emulate detection hardware
switch parameters.detection
    
    case 'quadrature'
        
        % Do nothing
        
    case 'uniaxial'
        
        % Destroy imaginary part
        fid=real(fid);
        
end

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(numel(parameters.sweep)~=1)||...
   (~isreal(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real number.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'drop_field')
    error('parameters.drop_field field is missing.');
elseif (~isnumeric(parameters.drop_field))||(~isreal(parameters.drop_field))||...
       (numel(parameters.drop_field)~=1)||(parameters.drop_field<0)
    error('parameters.drop_field should be a non-negative real scalar.');
end
if ~isfield(parameters,'drop_time')
    error('parameters.drop_time field is missing.');
elseif (~isnumeric(parameters.drop_time))||(~isreal(parameters.drop_time))||...
       (numel(parameters.drop_time)~=1)||(parameters.drop_time<=0)
    error('parameters.drop_time should be a positive real scalar.');
end
if ~isfield(parameters,'drop_npoints')
    error('parameters.drop_npoints field is missing.');
elseif (~isnumeric(parameters.drop_npoints))||(~isreal(parameters.drop_npoints))||...
       (numel(parameters.drop_npoints)~=1)||(parameters.drop_npoints<1)||...
       (mod(parameters.drop_npoints,1)~=0)
    error('parameters.drop_npoints should be a positive integer.');
end
if ~isfield(parameters,'drop_rate')
    error('parameters.drop_rate field is missing.');
elseif (~isnumeric(parameters.drop_rate))||(~isreal(parameters.drop_rate))||...
       (numel(parameters.drop_rate)~=1)||(parameters.drop_rate<=0)
    error('parameters.drop_rate should be a positive real scalar.');
end
if ~isfield(parameters,'detection')
    error('detection mode must be specified in parameters.detection');
end
if ~ismember(parameters.detection,{'uniaxial','quadrature'})
    error('parameters.detection must be ''uniaxial'' or ''quadrature''');
end
if ~isfield(parameters,'flip_angle')
    error('proton flip angle must be specified in parameters.flip_angle');
elseif (~isnumeric(parameters.flip_angle))||(~isreal(parameters.flip_angle))||...
       (numel(parameters.flip_angle)~=1)
    error('parameters.flip_angle should be a real scalar.');
end
end

% It wasn't a dark and stormy night. It should have been, but there's
% the weather for you. For every mad scientist who's had a convenient
% thunderstorm just on the night his Great Work is complete and lying
% on the slab, there have been dozens who've sat around aimlessly un-
% der the peaceful stars while Igor clocks up the overtime.
%
% Terry Pratchett


