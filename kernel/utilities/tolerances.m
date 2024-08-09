% Tolerances and fundamental constants. Sets various accuracy cut-offs,
% constants and tolerances used by Spinach kernel. Syntax:
%
%                 spin_system=tolerances(spin_system,sys)
%
% Parameters:
%
%    spin_system  -  Spinach system description object
%
%    sys          -  system specification object described
%                    in the input preparation section of 
%                    the manual
%
% Outputs:
%
%    spin_system  -  updated system description object
%
%    sys          -  system specification structure with
%                    the tolerance substructure parsed out
%
% Notes: direct calls and modifications to this function are discouraged:
%        the accuracy settings should be modified by setting the sys.tols
%        structure, see the input preparation manual.
%
% i.kuprov@soton.ac.uk
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=tolerances.m>

function [spin_system,sys]=tolerances(spin_system,sys)

% Check consistency
grumble(spin_system,sys);

% Interaction tensor clean-up tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')
    spin_system.tols.inter_cutoff=sys.tols.inter_cutoff; sys.tols=rmfield(sys.tols,'inter_cutoff');
    report(spin_system,[pad('Drop interaction tensors with 2-norms below (Hz)',65) ...
                        pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.inter_cutoff=eps();
    report(spin_system,[pad('Drop interaction tensors with 2-norms below (Hz)',65) ...
                        pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.inter_cutoff=1e-2;
    report(spin_system,[pad('Drop interaction tensors with 2-norms below (Hz)',65) ...
                        pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.inter_cutoff=1e-10;
    report(spin_system,[pad('Drop interaction tensors with 2-norms below (Hz)',65) ...
                        pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (safe default)']);
end

% Liouvillian matrix element zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')
    spin_system.tols.liouv_zero=sys.tols.liouv_zero; sys.tols=rmfield(sys.tols,'liouv_zero');
    report(spin_system,[pad('Drop Liouvillian elements with norm below (rad/s)',65) ...
                        pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.liouv_zero=eps();
    report(spin_system,[pad('Drop Liouvillian elements with norm below (rad/s)',65) ...
                        pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.liouv_zero=1e-5;
    report(spin_system,[pad('Drop Liouvillian elements with norm below (rad/s)',65) ...
                        pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.liouv_zero=1e-10;
    report(spin_system,[pad('Drop Liouvillian elements with norm below (rad/s)',65) ...
                        pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (safe default)']);
end

% Relaxation superoperator element zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_zero')
    spin_system.tols.rlx_zero=sys.tols.rlx_zero; sys.tols=rmfield(sys.tols,'rlx_zero');
    report(spin_system,[pad('Drop relaxation superoperator elements with norm below (Hz)',65) ...
                        pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.rlx_zero=eps();
    report(spin_system,[pad('Drop relaxation superoperator elements with norm below (Hz)',65) ...
                        pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.rlx_zero=1e-5;
    report(spin_system,[pad('Drop relaxation superoperator elements with norm below (Hz)',65) ...
                        pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.rlx_zero=1e-10;
    report(spin_system,[pad('Drop relaxation superoperator elements with norm below (Hz)',65) ...
                        pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (safe default)']);
end

% Propagator matrix element zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')
    spin_system.tols.prop_chop=sys.tols.prop_chop; sys.tols=rmfield(sys.tols,'prop_chop');
    report(spin_system,[pad('Drop propagator elements with norm below',65) ...
                        pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.prop_chop=eps();
    report(spin_system,[pad('Drop propagator elements with norm below',65) ...
                        pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.prop_chop=1e-8;
    report(spin_system,[pad('Drop propagator elements with norm below',65) ...
                        pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.prop_chop=1e-10;
    report(spin_system,[pad('Drop propagator elements with norm below',65) ...
                        pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (safe default)']);
end

% Subspace population tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'subs_drop')
    spin_system.tols.subs_drop=sys.tols.subs_drop; sys.tols=rmfield(sys.tols,'subs_drop');
    report(spin_system,[pad('Drop subspaces with population below',65) ...
                        pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.subs_drop=eps();
    report(spin_system,[pad('Drop subspaces with populations below',65) ...
                        pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.subs_drop=1e-2;
    report(spin_system,[pad('Drop subspaces with populations below',65) ...
                        pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.subs_drop=1e-10;
    report(spin_system,[pad('Drop subspaces with populations below',65) ...
                        pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (safe default)']);
end

% Irrep population tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')
    spin_system.tols.irrep_drop=sys.tols.irrep_drop; sys.tols=rmfield(sys.tols,'irrep_drop');
    report(spin_system,[pad('Drop irreps with populations below',65) ...
                        pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.irrep_drop=eps();
    report(spin_system,[pad('Drop irreps with populations below',65) ...
                        pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.irrep_drop=1e-2;
    report(spin_system,[pad('Drop irreps with populations below',65) ...
                        pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.irrep_drop=1e-10;
    report(spin_system,[pad('Drop irreps with populations below',65) ...
                        pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (safe default)']);
end

% ZTE sample length
if isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')
    spin_system.tols.zte_nsteps=sys.tols.zte_nsteps; sys.tols=rmfield(sys.tols,'zte_nsteps');
    report(spin_system,[pad('Number of time steps in ZTE sample',65) ...
                        pad(num2str(spin_system.tols.zte_nsteps),20) ' (user-specified)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.zte_nsteps=16;
    report(spin_system,[pad('Number of time steps in ZTE sample',65) ...
                        pad(num2str(spin_system.tols.zte_nsteps),20) ' (loose)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.zte_nsteps=NaN;
    report(spin_system,[pad('Number of time steps in ZTE sample',65) ...
                        pad('ZTE disabled',20) ' (paranoid)']);
else
    spin_system.tols.zte_nsteps=32;
    report(spin_system,[pad('Number of time steps in ZTE sample',65) ...
                        pad(num2str(spin_system.tols.zte_nsteps),20) ' (safe default)']);
end

% ZTE zero track tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')
    spin_system.tols.zte_tol=sys.tols.zte_tol; sys.tols=rmfield(sys.tols,'zte_tol');
    report(spin_system,[pad('ZTE track drop threshold',65) ...
                        pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (user-specified)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.zte_tol=1e-6;
    report(spin_system,[pad('ZTE track drop threshold',65) ...
                        pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (loose)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.zte_tol=NaN;
    report(spin_system,[pad('ZTE track drop threshold',65) ...
                        pad('ZTE disabled',20) ' (paranoid)']);
else
    spin_system.tols.zte_tol=1e-24;
    report(spin_system,[pad('ZTE track drop threshold',65) ...
                        pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (safe default)']);
end

% ZTE state vector density threshold 
if isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')
    spin_system.tols.zte_maxden=sys.tols.zte_maxden; sys.tols=rmfield(sys.tols,'zte_maxden');
    report(spin_system,[pad('ZTE off for state vector densities above',65) ...
                        pad(num2str(spin_system.tols.zte_maxden),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.zte_maxden=NaN;
    report(spin_system,[pad('ZTE off for state vector densities above',65) ...
                        pad('ZTE disabled',20) ' (paranoid)']);
else
    spin_system.tols.zte_maxden=0.5;
    report(spin_system,[pad('ZTE off for state vector densities above',65) ...
                        pad(num2str(spin_system.tols.zte_maxden),20) ' (safe default)']);
end

% Proximity tolerance for dipolar couplings (TODO: replace with energy tolerance)
if isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')
    spin_system.tols.prox_cutoff=sys.tols.prox_cutoff; sys.tols=rmfield(sys.tols,'prox_cutoff');
    report(spin_system,[pad('Drop dipolar couplings over more than (Angstrom)',65) ...
                        pad(num2str(spin_system.tols.prox_cutoff),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.prox_cutoff=inf();
    report(spin_system,[pad('Drop dipolar couplings over more than (Angstrom)',65) ...
                        pad(num2str(spin_system.tols.prox_cutoff),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.prox_cutoff=3.5;
    report(spin_system,[pad('Drop dipolar couplings over more than (Angstrom)',65) ...
                        pad(num2str(spin_system.tols.prox_cutoff),20) ' (loose)']);
else
    spin_system.tols.prox_cutoff=100;
    report(spin_system,[pad('Drop dipolar couplings over more than (Angstrom)',65) ...
                        pad(num2str(spin_system.tols.prox_cutoff),20) ' (safe default)']);
end

% Krylov method switchover
if isfield(sys,'tols')&&isfield(sys.tols,'krylov_tol')
    spin_system.tols.krylov_tol=sys.tols.krylov_tol; sys.tols=rmfield(sys.tols,'krylov_tol');
    report(spin_system,[pad('Krylov propagation for dimensions above',65) ...
                        pad(num2str(spin_system.tols.krylov_tol),20) ' (user-specified)']);    
else
    spin_system.tols.krylov_tol=10000;
    report(spin_system,[pad('Krylov propagation for dimensions above',65) ...
                        pad(num2str(spin_system.tols.krylov_tol),20) ' (safe default)']);
end

% Basis printing hush tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')
    spin_system.tols.basis_hush=sys.tols.basis_hush; sys.tols=rmfield(sys.tols,'basis_hush');
    report(spin_system,[pad('Basis set printing hush threshold',65) ...
                        pad(num2str(spin_system.tols.basis_hush),20) ' (user-specified)']);
else
    spin_system.tols.basis_hush=256;
    report(spin_system,[pad('Basis set printing hush threshold',65) ...
                        pad(num2str(spin_system.tols.basis_hush),20) ' (safe default)']);
end

% Subspace bundle size
if isfield(sys,'tols')&&isfield(sys.tols,'merge_dim')
    spin_system.tols.merge_dim=sys.tols.merge_dim; sys.tols=rmfield(sys.tols,'merge_dim');
    report(spin_system,[pad('Collect small subspaces into bundles of dimension',65) ...
                        pad(num2str(spin_system.tols.merge_dim),20) ' (user-specified)']);
else
    spin_system.tols.merge_dim=1000;
    report(spin_system,[pad('Collect small subspaces into bundles of dimension',65) ...
                        pad(num2str(spin_system.tols.merge_dim),20) ' (safe default)']);
end

% Sparse algebra tolerance on density
if isfield(sys,'tols')&&isfield(sys.tols,'dense_matrix')
    spin_system.tols.dense_matrix=sys.tols.dense_matrix; sys.tols=rmfield(sys.tols,'dense_matrix');
    report(spin_system,[pad('Force sparse algebra for matrix density below',65) ...
                        pad(num2str(spin_system.tols.dense_matrix),20) ' (user-specified)']);
else
    spin_system.tols.dense_matrix=0.15;
    report(spin_system,[pad('Force sparse algebra for matrix density below',65) ...
                        pad(num2str(spin_system.tols.dense_matrix),20) ' (safe default)']);
end

% Sparse algebra tolerance on dimension
if isfield(sys,'tols')&&isfield(sys.tols,'small_matrix')
    spin_system.tols.small_matrix=sys.tols.small_matrix; sys.tols=rmfield(sys.tols,'small_matrix');
    report(spin_system,[pad('Force full algebra for matrix dimension below',65) ...
                        pad(num2str(spin_system.tols.small_matrix),20) ' (user-specified)']);
else
    spin_system.tols.small_matrix=200;
    report(spin_system,[pad('Force full algebra for matrix dimension below',65) ...
                        pad(num2str(spin_system.tols.small_matrix),20) ' (safe default)']);
end

% Relative accuracy of the elements of Redfield superoperator
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_integration')
    spin_system.tols.rlx_integration=sys.tols.rlx_integration; sys.tols=rmfield(sys.tols,'rlx_integration');
    report(spin_system,[pad('Rel. acc. for the elements of Redfield superoperator',65) ...
                        pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (user-specified)']);
elseif ismember('paranoia',spin_system.sys.enable)
    spin_system.tols.rlx_integration=1e-6;
    report(spin_system,[pad('Rel. acc. for the elements of Redfield superoperator',65) ...
                        pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (paranoid)']);
elseif ismember('cowboy',spin_system.sys.enable)
    spin_system.tols.rlx_integration=1e-2;
    report(spin_system,[pad('Rel. acc. for the elements of Redfield superoperator',65) ...
                        pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.rlx_integration=1e-4;
    report(spin_system,[pad('Rel. acc. for the elements of Redfield superoperator',65) ...
                        pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (safe default)']);
end

% Algorithm selection for propagator derivatives
if isfield(sys,'tols')&&isfield(sys.tols,'dP_method')
    spin_system.tols.dP_method=sys.tols.dP_method; sys.tols=rmfield(sys.tols,'dP_method');
    report(spin_system,[pad('Matrix exponential differentiation algorithm',65) ...
                        pad(spin_system.tols.dP_method,20) ' (user-specified)']);
else
    spin_system.tols.dP_method='auxmat';
    report(spin_system,[pad('Matrix exponential differentiation algorithm',65) ...
                        pad(spin_system.tols.dP_method,20) ' (safe default)']);
end

% Number of PBC images for dipolar couplings
if isfield(sys,'tols')&&isfield(sys.tols,'dd_ncells')
    spin_system.tols.dd_ncells=sys.tols.dd_ncells; sys.tols=rmfield(sys.tols,'dd_ncells');
    report(spin_system,[pad('Number of PBC images for DD couplings in PBC systems',65) ...
                        pad(num2str(spin_system.tols.dd_ncells),20) ' (user-specified)']);
else
    spin_system.tols.dd_ncells=2;
    report(spin_system,[pad('Number of PBC images for DD couplings in PBC systems',65) ...
                        pad(num2str(spin_system.tols.dd_ncells),20) ' (safe default)']);
end

% Cache storage timeout
if isfield(sys,'tols')&&isfield(sys.tols,'cache_mem')
    spin_system.tols.cache_mem=sys.tols.cache_mem; sys.tols=rmfield(sys.tols,'cache_mem');
    report(spin_system,[pad('Number of days before a cache record is deleted',65) ...
                        pad(num2str(spin_system.tols.cache_mem),20) ' (user-specified)']);
else
    spin_system.tols.cache_mem=30;
    report(spin_system,[pad('Number of days before a cache record is deleted',65) ...
                        pad(num2str(spin_system.tols.cache_mem),20) ' (safe default)']);
end

% Fundamental constants
spin_system.tols.hbar=6.62607015e-34/(2*pi); % J*s, exact number
report(spin_system,[pad('Planck constant (hbar)',65) pad(num2str(spin_system.tols.hbar,'%0.8e'),20)]);
spin_system.tols.kbol=1.380649e-23; % J*K^{-1}, exact number 
report(spin_system,[pad('Boltzmann constant (k)',65) pad(num2str(spin_system.tols.kbol,'%0.8e'),20)]);
spin_system.tols.freeg=2.00231930436256;
report(spin_system,[pad('Free electron g-factor',65) pad(num2str(spin_system.tols.freeg,'%0.8e'),20)]);
spin_system.tols.mu0=1.25663706212e-6; % N A^-2 = T^2 m^3 J^-1
report(spin_system,[pad('Vacuum permeability',65) pad(num2str(spin_system.tols.mu0,'%0.8e'),20)]);
spin_system.tols.muB=9.2740100783e-24; % J T^-1
report(spin_system,[pad('Bohr magneton',65) pad(num2str(spin_system.tols.muB,'%0.8e'),20)]);

% Paranoia switches
if ismember('paranoia',spin_system.sys.enable)

    % Make sure zero track elimination is disabled
    spin_system.sys.disable=unique([spin_system.sys.disable {'zte'}]);

    % Make sure operator and propagator caching is not enabled
    spin_system.sys.enable=setdiff(spin_system.sys.enable,{'op_cache','prop_cache'});
    
end

% Catch unparsed options
if isfield(sys,'tols')
    unparsed=fieldnames(sys.tols);
    if ~isempty(unparsed)
        for n=1:numel(unparsed)
            report(spin_system,['unrecognised option - ' unparsed{n}]);
        end
        error('unrecognised options in sys.tols');
    end
    sys=rmfield(sys,'tols');
end

end

% Consistency enforcement
function grumble(spin_system,sys) %#ok<INUSL>
if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')
    if (~isnumeric(sys.tols.inter_cutoff))||(~isreal(sys.tols.inter_cutoff))||...
       (~isscalar(sys.tols.inter_cutoff))||(sys.tols.inter_cutoff<0)
        error('sys.tols.inter_cutoff must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')
    if (~isnumeric(sys.tols.liouv_zero))||(~isreal(sys.tols.liouv_zero))||...
       (~isscalar(sys.tols.liouv_zero))||(sys.tols.liouv_zero<0)
        error('sys.tols.liouv_zero must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')
    if (~isnumeric(sys.tols.prop_chop))||(~isreal(sys.tols.prop_chop))||...
       (~isscalar(sys.tols.prop_chop))||(sys.tols.prop_chop<0)
        error('sys.tols.prop_chop must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'subs_drop')
    if (~isnumeric(sys.tols.subs_drop))||(~isreal(sys.tols.subs_drop))||...
       (~isscalar(sys.tols.subs_drop))||(sys.tols.subs_drop<0)
        error('sys.tols.subs_drop must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')
    if (~isnumeric(sys.tols.irrep_drop))||(~isreal(sys.tols.irrep_drop))||...
       (~isscalar(sys.tols.irrep_drop))||(sys.tols.irrep_drop<0)
        error('sys.tols.irrep_drop must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'path_drop')
    if (~isnumeric(sys.tols.path_drop))||(~isreal(sys.tols.path_drop))||...
       (~isscalar(sys.tols.path_drop))||(sys.tols.path_drop<0)
        error('sys.tols.path_drop must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')
    if (~isnumeric(sys.tols.zte_tol))||(~isreal(sys.tols.zte_tol))||...
       (~isscalar(sys.tols.zte_tol))||(sys.tols.zte_tol<0)
        error('sys.tols.zte_tol must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')
    if (~isnumeric(sys.tols.zte_nsteps))||(~isreal(sys.tols.zte_nsteps))||...
       (~isscalar(sys.tols.zte_nsteps))||(mod(sys.tols.zte_nsteps,1)~=0)||...
       (sys.tols.zte_nsteps<1)
        error('sys.tols.zte_tol must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')
    if (~isnumeric(sys.tols.zte_maxden))||(~isreal(sys.tols.zte_maxden))||...
       (~isscalar(sys.tols.zte_maxden))||(sys.tols.zte_maxden<0)||...
       (sys.tols.zte_maxden>1)
        error('sys.tols.zte_maxden must be a real number between 0 and 1.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'dense_matrix')
    if (~isnumeric(sys.tols.dense_matrix))||(~isreal(sys.tols.dense_matrix))||...
       (~isscalar(sys.tols.dense_matrix))||(sys.tols.dense_matrix<0)||...
       (sys.tols.dense_matrix>1)
        error('sys.tols.dense_matrix must be a real number between 0 and 1.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_integration')
    if (~isnumeric(sys.tols.rlx_integration))||(~isreal(sys.tols.rlx_integration))||...
       (~isscalar(sys.tols.rlx_integration))||(sys.tols.rlx_integration<=0)
        error('sys.tols.rlx_integration must be a positive real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')
    if (~isnumeric(sys.tols.prox_cutoff))||(~isreal(sys.tols.prox_cutoff))||...
       (~isscalar(sys.tols.prox_cutoff))||(sys.tols.prox_cutoff<0)
        error('sys.tols.prox_cutoff must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'krylov_tol')
    if (~isnumeric(sys.tols.krylov_tol))||(~isreal(sys.tols.krylov_tol))||...
       (~isscalar(sys.tols.krylov_tol))||(mod(sys.tols.krylov_tol,1)~=0)||...
       (sys.tols.krylov_tol<1)
        error('sys.tols.krylov_tol must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')
    if (~isnumeric(sys.tols.basis_hush))||(~isreal(sys.tols.basis_hush))||...
       (~isscalar(sys.tols.basis_hush))||(mod(sys.tols.basis_hush,1)~=0)||...
       (sys.tols.basis_hush<1)
        error('sys.tols.basis_hush must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'merge_dim')
    if (~isnumeric(sys.tols.merge_dim))||(~isreal(sys.tols.merge_dim))||...
       (~isscalar(sys.tols.merge_dim))||(mod(sys.tols.merge_dim,1)~=0)||...
       (sys.tols.merge_dim<1)
        error('sys.tols.merge_dim must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'small_matrix')
    if (~isnumeric(sys.tols.small_matrix))||(~isreal(sys.tols.small_matrix))||...
       (~isscalar(sys.tols.small_matrix))||(mod(sys.tols.small_matrix,1)~=0)||...
       (sys.tols.small_matrix<1)
        error('sys.tols.small_matrix must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'dd_ncells')
    if (~isnumeric(sys.tols.dd_ncells))||(~isreal(sys.tols.dd_ncells))||...
       (~isscalar(sys.tols.dd_ncells))||(mod(sys.tols.dd_ncells,1)~=0)||...
       (sys.tols.dd_ncells<1)
        error('sys.tols.dd_ncells must be a positive real integer.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'cache_mem')
    if (~isnumeric(sys.tols.cache_mem))||(~isreal(sys.tols.cache_mem))||...
       (~isscalar(sys.tols.cache_mem))||(sys.tols.cache_mem<0)
        error('sys.tols.cache_mem must be a non-negative real number.');
    end
end
if isfield(sys,'tols')&&isfield(sys.tols,'dP_method')
    if ~ischar(sys.tols.dP_method)
        error('sys.tols.dP_method must be a character string.');
    end
end
end

% Man once surrendering his reason, has no remaining guard against
% absurdities the most monstrous, and like a ship without rudder, 
% is the sport of every wind. 
%
% Thomas Jefferson

