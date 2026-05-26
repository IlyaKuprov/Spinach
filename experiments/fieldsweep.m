% Field swept powder EPR spectra. A rough implementation with ex-
% pensive eigenfields algorithm, an explicit spherical grid, and
% a hard-coded Lorentzian line shape. Syntax:
%
%        [b_axis,spec]=fieldsweep(spin_system,parameters)
%
% Parameters:
%
%     parameters.grid      -  initial spherical grid, ideally
%                             a non-symmetric one to avoid
%                             transition degeneracies, a good
%                             start is
%                             'rep_2ang_100pts_sph'
%
%     parameters.spins     -  a one-element cell array speci-
%                             fying the spin that is coupled
%                             couple to the microwave field,
%                             e.g. {'E'}
%
%     parameters.mw_freq   -  microwave frequency, Hz
%
%     parameters.fwhm      -  Lorentzian line FWHM, Tesla
%
%     parameters.window    -  field sweep window in Tesla,
%                             as a vector [Bmin Bmax]
%
%     parameters.npoints   -  number of points in the sweep
%
%     parameters.tm_tol    -  relative transition moment to-
%                             lerance, 0.01 is a good start
%
%     parameters.rspt_order - perturbation theory order for
%                             eigenfields calculation, 2 is
%                             a good start; specify Inf for
%                             exact diagonalisation
%
%     parameters.int_tol   -  powder integration tolerance,
%                             a balance between speed and
%                             integration accuracy
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

function [spec,parameters]=fieldsweep(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Peak position tolerance is quarter a grid pixel
window_size=max(parameters.window)-min(parameters.window);
parameters.pp_tol=0.25*window_size/(parameters.npoints-1);

% Get the Hamiltonians and tidy up their isotropic parts
[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));
[Iz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));
Ic=(Ic+Ic')/2; Iz=(Iz+Iz')/2; 

% Get the microwave operator (state in Liouville space)
Hmw=state(spin_system,'Lx',parameters.spins{1});

% Strip the spin system object down to minimum size
ss_parfor=stripper(spin_system,'field_sweep');

% Get the initial grid and compute its convex hull
load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' ... 
      filesep parameters.grid],'betas','gammas');
hull=get_hull(betas,gammas); grid_size=numel(betas);

% Make the magnetic field axis
parameters.b_axis=linspace(parameters.window(1),...
                           parameters.window(2),...
                           parameters.npoints);

% Preallocate asynchronous work outputs
eigensets(1:grid_size)=parallel.Future;

% Schedule work
for n=1:grid_size

    % Create a local copy and specify system orientation 
    localpar=parameters; localpar.orientation=[0 betas(n) gammas(n)];

    % Assemble Zeeman and coupling Hamiltonians
    Hz=Iz+orientation(Qz,localpar.orientation); 
    Hc=Ic+orientation(Qc,localpar.orientation); 
    
    % Submit asynchronous eigenfields calculation job for this vertex
    eigensets(n)=parfeval(@eigenfields,1,ss_parfor,localpar,Hz,Hc,Hmw);

end

% Retrieve vertex eigensets
eigensets=fetchOutputs(eigensets);

% Add coordinates
for n=1:grid_size
    eigensets(n).xyz=[sin(betas(n))*cos(gammas(n)); 
                      sin(betas(n))*sin(gammas(n)); 
                      cos(betas(n))];
end

% Preallocate asynchronous work outputs
spec(1:size(hull,1))=parallel.Future;

% Schedule work 
for n=1:size(hull,1)
    
    % Extract triangle vertex indices
    a=hull(n,1); b=hull(n,2); c=hull(n,3);

    % Extract vertex eigensets
    triangle=eigensets([a b c]); 
    
    % Call recursive Voitlander integrator
    spec(n)=parfeval(@voitlander,1,ss_parfor,parameters,...
                                   triangle,Ic,Iz,Qc,Qz,Hmw);

end

% Retrieve and sum up cell integrals
spec=sum(fetchOutputs(spec),1);

end

% Consistency enforcement
function grumble(spin_system,parameters)
if spin_system.inter.magnet~=1
    error('field swept experiment, set sys.magnet=1');
end
if ~isfield(parameters,'grid')
    error('grid must be specified in parameters.grid variable.');
end
if ~ischar(parameters.grid)
    error('parameters.grid must be character string.');
end
if ~isfield(parameters,'spins')
    error('spin to couple to the microwave field must be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
   (~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
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
   (~isreal(parameters.npoints))||(parameters.npoints<=1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a real integer greater than one.');
end
if ~isfield(parameters,'window')
    error('sweep window should be specified in parameters.window variable.');
end
if (~isnumeric(parameters.window))||(numel(parameters.window)~=2)||...
   (~isreal(parameters.window))||(parameters.window(2)<=parameters.window(1))
    error('parameters.window must contain two ascending numbers.');
end
if ~isfield(parameters,'tm_tol')
    error('transition moment tolerance must be specified in parameters.tm_tol variable.');
end
if (~isnumeric(parameters.tm_tol))||(~isscalar(parameters.tm_tol))||...
   (~isreal(parameters.tm_tol))||(parameters.tm_tol<0)
    error('parameters.tm_tol must be a non-negative real scalar.');
end
if ~isfield(parameters,'rspt_order')
    error('perturbation theory order must be specified in parameters.rspt_order variable.');
end
if (~isnumeric(parameters.rspt_order))||(~isscalar(parameters.rspt_order))||...
   (~isreal(parameters.rspt_order))||((mod(parameters.rspt_order,1)~=0)&&...
   (~isinf(parameters.rspt_order)))||(parameters.rspt_order<0)
    error('parameters.rspt_order must be a non-negative integer or Inf.');
end
if ~isfield(parameters,'int_tol')
    error('integration tolerance must be specified in parameters.int_tol variable.');
end
if (~isnumeric(parameters.int_tol))||(~isscalar(parameters.int_tol))||...
   (~isreal(parameters.int_tol))||(parameters.int_tol<=0)
    error('parameters.int_tol must be a positive real scalar.');
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
