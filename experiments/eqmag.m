% Computes the molar magnetization vector at the thermal equilibrium at 
% the temperature specified in inter.temperature and magnetic field spe-
% cified in sys.magnet (assumed to be along the Z-axis), averaged over 
% system orientations using the spherical grid specified. Syntax:
%
%                   magn=eqmag(spin_system,parameters)
%
% Input:
%
%       parameters.grid - spherical grid for averaging
%
% Output:
%
%       magn - molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
%
% Note: the use of bas.formalism='zeeman-hilb' is obligatory.
%
% Note: Spinach uses NMR convention for the exchange coupling: exchange
%       interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
%       where J is in Hz.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=eqmag.m>

function magn=eqmag(spin_system,parameters)

% Check consistency
grumble(parameters)

% Get the number of spins
nspins=spin_system.comp.nspins;

% Get the g-tensor for each spin in Bohr magneton units
for n=nspins:-1:1
    g{n}=gtensorof(spin_system,n);  
end

% Get Sx Sy Sz operators for each spin
for n=nspins:-1:1
    Sx{n}=(operator(spin_system,{'L+'},{n})+operator(spin_system,{'L-'},{n}))/2;
    Sy{n}=(operator(spin_system,{'L+'},{n})-operator(spin_system,{'L-'},{n}))/2i;
    Sz{n}=operator(spin_system,{'Lz'},{n});
end

% Define the problem dimension
problem_dim=prod(spin_system.comp.mults);

% Get Hamiltonian components
[I,Q]=hamiltonian(assume(spin_system,'labframe'));

% Load the spherical grid
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' parameters.grid '.mat']);
weights=sph_grid.weights; alphas=sph_grid.alphas;
betas=sph_grid.betas; gammas=sph_grid.gammas;

% Get the magnetization vector started
magn=[0 0 0];

% Loop over the grid
parfor k=1:numel(weights) %#ok<*PFBNS>
    
    % Get the rotation matrix
    R=euler2dcm([alphas(k) betas(k) gammas(k)]);
    
    % Start magnetic moment operators
    mx=sparse(problem_dim,problem_dim);
    my=sparse(problem_dim,problem_dim);
    mz=sparse(problem_dim,problem_dim);
    
    % Compute the detection operator
    for n=1:nspins
        
        % Rotate g-tensor (because we rotate the molecule, not the field)
        g_rot=sparse(R*g{n}*R');
        
        % Build magnetic moment operators
        mx=mx-g_rot(1,1)*Sx{n}-g_rot(1,2)*Sy{n}-g_rot(1,3)*Sz{n}; 
        my=my-g_rot(2,1)*Sx{n}-g_rot(2,2)*Sy{n}-g_rot(2,3)*Sz{n};
        mz=mz-g_rot(3,1)*Sx{n}-g_rot(3,2)*Sy{n}-g_rot(3,3)*Sz{n}; 
        
    end
    
    % Get the equilibrium density matrix
    rho=equilibrium(spin_system,I,Q,[alphas(k) betas(k) gammas(k)]);
    
    % Normalize density
    rho=rho./trace(rho);

    % Get the magnetisation vector in Bohr magnetons per mol
    magn=magn+sph_grid.weights(k)*[trace(rho*mx) trace(rho*my) trace(rho*mz)];

end

% Keep the real part
magn=real(magn);

end

% Consistency checking
function grumble(parameters)
if ~isfield(parameters,'grid')
    error('spherical averaging grid must be specified in parameters.grid variable.');
elseif isempty(parameters.grid)
    error('parameters.grid variable cannot be empty.');
elseif ~ischar(parameters.grid)
    error('parameters.grid variable must be a character string.');
end
end

% The window to the world can be covered by a newspaper. 
%
% Stanislaw Jerzy Lec

