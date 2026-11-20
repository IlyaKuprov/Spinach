% Drift Hamiltonians under magic angle spinning, resolved in rotor
% phase, for every combination of a two-angle powder grid orienta-
% tion and an initial rotor phase. The Hamiltonian is taken to se-
% cond order in the rotating frame of the nucleus, and is held 
% constant within each rotor phase tick. Syntax:
%
%              [drifts,weights]=mas_drifts(spin_system,parameters)
%
% Parameters:
%
%    parameters.spins    - the nucleus, a cell array with one
%                          isotope string, e.g. {'27Al'}
%
%    parameters.rframes  - optional ordered {isotope,order} pairs;
%                          replaces the legacy spins specification
%
%    parameters.axis     - spinning axis, a normalised row
%                          vector with three elements
%
%    parameters.grid     - powder grid name or grid_check structure
%
%    parameters.n_ticks  - rotor phase ticks per rotor period
%
%    parameters.n_phases - number of initial rotor phases, must
%                          be a divisor of parameters.n_ticks
%
%    parameters.n_slices - number of ticks in the pulse
%
% Outputs:
%
%    drifts   - cell array over the ensemble, grid orientations in
%               the outer index and initial rotor phases in the in-
%               ner index; each element is a cell array of n_slices
%               Hamiltonian matrices, one per tick of the pulse
%
%    weights  - column of grid weights divided over initial phases
%
% Note: the crystallite orientation uses all three Euler angles of the
%       grid, as in singlerot.m; two-angle grids keep the azimuth in
%       the third angle and have zero first angles, which the rotor
%       phase then supplies. The spin system must carry laboratory
%       frame assumptions, call assume() first.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mas_drifts.m>

function [drifts,weights]=mas_drifts(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Get the laboratory frame Hamiltonian
[H,Q]=hamiltonian(spin_system);

% Resolve the legacy single frame or the ordered numerical frames
if isfield(parameters,'rframes')
    frames=parameters.rframes;
else
    frames={{parameters.spins{1},2}};
end
C=cellfun(@(frame)carrier(spin_system,frame{1}),frames,'UniformOutput',false);

% Remove all active carriers before reintroducing each as the reference
if numel(frames)>1
    for f=1:numel(frames)
        H=H-C{f};
    end
end

% Read the integration grid and preserve its quadrature weights
if ischar(parameters.grid)
    sph_grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',parameters.grid),...
                  'alphas','betas','gammas','weights');
else
    sph_grid=parameters.grid;
end
grid_check(sph_grid);
if any(sph_grid.weights<0)||~any(sph_grid.weights>0)
    error('the powder grid must have nonnegative weights and positive total weight.');
end
if (nargout<2)&&any(sph_grid.weights~=sph_grid.weights(1))
    error('Spinach:mas_drifts:weightsRequired','request the weights output for nonuniform powder grids.');
end
weights=repelem(sph_grid.weights(:)/parameters.n_phases,parameters.n_phases);

% Get rotor axis orientation
[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1),...
                                   parameters.axis(2),...
                                   parameters.axis(3));
rotor_theta=pi/2-rotor_theta;

% Rotor phases at tick midpoints
rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks;

% Tick shift between initial rotor phases
tick_shift=parameters.n_ticks/parameters.n_phases;

% Silence the workers
spin_system.sys.output='hush';

% Parallel loop over grid orientations
n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients);
parfor n=1:n_orients %#ok<*PFBNS>

    % Preallocate the rotor stack
    stack=cell(1,parameters.n_ticks);

    % Loop over rotor phase ticks
    for k=1:parameters.n_ticks

        % Start with the isotropic part
        stack{k}=H;

        % Loop over spherical ranks
        for r=1:numel(Q)

            % Compute crystallite orientation
            D_mol2rot=wigner(r,sph_grid.alphas(n),sph_grid.betas(n),sph_grid.gammas(n));

            % Compute rotor axis tilt
            D_lab2rot=wigner(r,rotor_phi,rotor_theta,0);

            % Compute rotor rotation
            D_rotor=wigner(r,0,0,rotor_phases(k));

            % Compose rotations
            D_comp=D_lab2rot*D_rotor*D_mol2rot;

            % Build the anisotropic part
            for p=1:(2*r+1)
                for q=1:(2*r+1)
                    stack{k}=stack{k}+D_comp(p,q)*Q{r}{p,q};
                end
            end

        end

        % Apply numerical frames in the supplied order
        stack{k}=(stack{k}+stack{k}')/2;
        for f=1:numel(frames)
            if numel(frames)>1
                stack{k}=stack{k}+C{f};
            end
            stack{k}=rotframe(spin_system,C{f},(stack{k}+stack{k}')/2,...
                              frames{f}{1},frames{f}{2});
        end
        stack{k}=full(stack{k});

    end

    % Slice the rotor stack for each initial rotor phase
    members=cell(1,parameters.n_phases);
    for j=1:parameters.n_phases
        tick_idx=mod((0:(parameters.n_slices-1))+(j-1)*tick_shift,...
                     parameters.n_ticks)+1;
        members{j}=stack(tick_idx);
    end
    orient_drifts{n}=members;

end

% Concatenate the ensemble
drifts=[orient_drifts{:}];

end

% Consistency enforcement
function grumble(spin_system,parameters)
if isfield(parameters,'rframes')
    frames=parameters.rframes;
    if ~iscell(frames)||isempty(frames)||~isvector(frames)
        error('parameters.rframes must be a nonempty cell vector.');
    end
    isotopes=cell(size(frames));
    for n=1:numel(frames)
        frame=frames{n};
        if ~iscell(frame)||numel(frame)~=2||~ischar(frame{1})||...
           ~isnumeric(frame{2})||~isreal(frame{2})||~isscalar(frame{2})||...
           isnan(frame{2})||frame{2}<0||(isfinite(frame{2})&&mod(frame{2},1)~=0)
            error('each frame must be an isotope and a nonnegative integer or Inf order.');
        end
        isotopes{n}=frame{1};
    end
    if numel(unique(isotopes))~=numel(isotopes)||~all(ismember(isotopes,spin_system.comp.isotopes))
        error('frame isotopes must be distinct and present in the spin system.');
    end
elseif ~isfield(parameters,'spins')||~iscell(parameters.spins)||...
       numel(parameters.spins)~=1||~ischar(parameters.spins{1})||...
       ~ismember(parameters.spins{1},spin_system.comp.isotopes)
    error('parameters.spins must contain one isotope present in the spin system.');
end
if ~isfield(parameters,'axis')
    error('spinning axis must be specified in parameters.axis field.');
end
if (~isnumeric(parameters.axis))||(~isreal(parameters.axis))||...
   (~isrow(parameters.axis))||(numel(parameters.axis)~=3)||...
   any(~isfinite(parameters.axis))||(abs(norm(parameters.axis,2)-1)>1e-6)
    error('parameters.axis must be a normalised row vector with three finite real elements.');
end
if ~isfield(parameters,'grid')
    error('powder grid must be specified in parameters.grid field.');
end
if ~ischar(parameters.grid)&&~isstruct(parameters.grid)
    error('parameters.grid must be a grid name or a grid structure.');
end
if ~isfield(parameters,'n_ticks')
    error('ticks per rotor period must be specified in parameters.n_ticks field.');
end
if (~isnumeric(parameters.n_ticks))||(~isreal(parameters.n_ticks))||...
   (~isscalar(parameters.n_ticks))||(parameters.n_ticks<1)||...
   (mod(parameters.n_ticks,1)~=0)
    error('parameters.n_ticks must be a positive integer.');
end
if ~isfield(parameters,'n_phases')
    error('number of initial rotor phases must be specified in parameters.n_phases field.');
end
if (~isnumeric(parameters.n_phases))||(~isreal(parameters.n_phases))||...
   (~isscalar(parameters.n_phases))||(parameters.n_phases<1)||...
   (mod(parameters.n_phases,1)~=0)||(mod(parameters.n_ticks,parameters.n_phases)~=0)
    error('parameters.n_phases must be a positive integer divisor of parameters.n_ticks.');
end
if ~isfield(parameters,'n_slices')
    error('number of pulse slices must be specified in parameters.n_slices field.');
end
if (~isnumeric(parameters.n_slices))||(~isreal(parameters.n_slices))||...
   (~isscalar(parameters.n_slices))||(parameters.n_slices<1)||...
   (mod(parameters.n_slices,1)~=0)
    error('parameters.n_slices must be a positive integer.');
end
end

% Everything should be made as simple as possible,
% but not simpler.
%
% Albert Einstein

