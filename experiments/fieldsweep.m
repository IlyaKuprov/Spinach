% Field swept powder EPR spectra. A very rough implementation using
% expensive eigenfields algorithm, an explicit spherical grid, and
% a hard-coded Lorentzian line shape function in Fermi Golden Rule.
%
%          [b_axis,spec]=fieldsweep(spin_system,parameters)
%
% Parameters:
%
%     parameters.grid      -  starting rank of the adaptively 
%                             subdivided spherical quadrature
%                             grid, 6 is a good choice
%
%     parameters.mw_freq   -  microwave frequency, Hz
%
%     parameters.fwhm      -  line FWHM, Tesla
%
%     parameters.window    -  field sweep window, [Bmin Bmax]
%
%     parameters.npoints   -  number of points in the sweep
%
% Outputs:
%
%     b_axis   - magnetic field axis for plotting
%
%     spec     - field-swept EPR spectrum
%
% Note: irrespective of the actual sweep extents, the magnetic field
%       in sys.magnet should be set to 1 Tesla.
%
% Note: this experiment should be called directly without a context.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fieldsweep.m>

function [b_axis,spec]=fieldsweep(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Set peak position tolerance to half-pixel
parameters.pp_tol=0.5*(max(parameters.window)-...
                       min(parameters.window))/(parameters.npoints-1);

% Get the Hamiltonians
[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));
[Iz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));

% Tidy up Hamiltonians
Ic=(Ic+Ic')/2; Iz=(Iz+Iz')/2; 

% Get the microwave operator (vector in Liouville space)
Hmw=(state(spin_system,'L+',parameters.spins{1})+...
     state(spin_system,'L-',parameters.spins{1}))/2;

% Get the initial grid and compute its convex hull
[alps,bets,gams]=grid_trian('stoll',parameters.grid); 
hull=get_hull(bets,gams); grid_size=numel(alps);

% Make the field axis
b_axis=linspace(parameters.window(1),...
                parameters.window(2),...
                parameters.npoints);

% Preallocate the spectrum
spec=zeros(size(b_axis),'like',1i);

% Eigenfields at grid vertices
tf=cell(grid_size,1); tm=cell(grid_size,1); 
tw=cell(grid_size,1); pd=cell(grid_size,1);
parfor n=1:grid_size %#ok<*PFBNS>

    % Localise parameters array and set the orientation
    loc_params=parameters; loc_params.orientation=[alps(n) bets(n) gams(n)];
    
    % Transition fields and moments
    [tf{n},tm{n},tw{n},pd{n}]=eigenfields(spin_system,loc_params,Iz,Qz,Ic,Qc,Hmw); 
    
end

% Convert grid to Cartesian coodinates
x=sin(bets).*cos(gams); y=sin(bets).*sin(gams); z=cos(bets);

% Voitlander integrator at each triangle
parfor n=1:size(hull,1)
    
    % Extract triangle indices
    a=hull(n,1); b=hull(n,2); c=hull(n,3);

    % Get vertex coordinates
    xyz_a=[x(a); y(a); z(a)]; 
    xyz_b=[x(b); y(b); z(b)]; 
    xyz_c=[x(c); y(c); z(c)];
    
    % Call Voitlander integrator
    spec=spec+voitlander(spin_system,parameters,xyz_a,xyz_b,xyz_c,...
                                     tf{a},tf{b},tf{c},tm{a},tm{b},tm{c},...
                                     tw{a},tw{b},tw{c},pd{a},pd{b},pd{c},...
                                     Ic,Iz,Qc,Qz,Hmw,b_axis);
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters)
if spin_system.inter.magnet~=1
    error('field swept experiment, set sys.magnet=1');
end
if ~isfield(parameters,'grid')
    error('stoll grid rank must be specified in parameters.grid variable.');
end
if (~isnumeric(parameters.grid))||(~isscalar(parameters.grid))||...
   (~isreal(parameters.grid))||(parameters.grid<=0)
    error('parameters.grid must be a positive real integer.');
end
if ~isfield(parameters,'mw_freq')
    error('MW frequency should be specified in parameters.mw_freq variable.');
end
if (~isnumeric(parameters.mw_freq))||(~isscalar(parameters.mw_freq))||...
   (~isreal(parameters.mw_freq))
    error('parameters.mw_freq must be a real scalar.');
end
if ~isfield(parameters,'fwhm')
    error('line width should be specified in parameters.fwhm variable.');
end
if (~isnumeric(parameters.fwhm))||(~isscalar(parameters.fwhm))||...
   (~isreal(parameters.fwhm))||(parameters.fwhm<=0)
    error('parameters.fwhm must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isscalar(parameters.npoints))||...
   (~isreal(parameters.npoints))||(parameters.npoints<=0)
    error('parameters.npoints must be a positive real integer.');
end
if ~isfield(parameters,'window')
    error('sweep window should be specified in parameters.window variable.');
end
if (~isnumeric(parameters.window))||(numel(parameters.window)~=2)||...
   (~isreal(parameters.window))||(parameters.window(2)<=parameters.window(1))
    error('parameters.window must contain two ascending numbers.');
end
end

% He leads the existence of a real bohemian intellectual. Washing, 
% grooming, and changing his linen are things he does rarely, and
% he likes to get drunk. Though he is often idle for days on end,
% he will work day and night with tireless endurance when he has 
% a great deal of work to do. He has no fixed times for going to 
% sleep and waking up. He often stays up all night, and then lies
% down fully clothed on the sofa at midday and sleeps till evening,
% untroubled by the comings and goings of the whole world.
%
% From the Prussian secret police
% report on Karl Marx

