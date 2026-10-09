% Extended T1/T2 relaxation model returning the relaxation super-
% operators separately for the longitudinal and the transverse 
% states. Syntax:
%
%       [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)
%
% Parameters:
%
%    euler_angles  - three Euler angles (ZYZ active convention
%                    in radians) specifying system orientation
%                    relative to the input orientation; requi-
%                    red when R1 and/or R2 rates had been spe-
%                    cified as 3x3 tensor or a function handle,
%                    this argument has no effect for R1 and R2
%                    rates specified as scalars.
%
% Outputs:
%
%    R1Op - relaxation superoperator containing
%           all longitudinal relaxation terms
%
%    R2Op - relaxation superoperator containing
%           all transverse relaxation terms
%
% Note: multi-spin orders relax at the sum of the rates of
%       their constituent single-spin orders.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_t1_t2.m>

function [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)

% Check consistency
grumble(spin_system);

% Preallocate relaxation rates
r1_rates=zeros(spin_system.comp.nspins,1);
r2_rates=zeros(spin_system.comp.nspins,1);

% Fill in R1 relaxation rates
for n=1:numel(spin_system.comp.isotopes)

    % Get the specification
    current_r1_rate=spin_system.rlx.r1_rates{n};

    % Scalar relaxation rate
    if isnumeric(current_r1_rate)&&isscalar(current_r1_rate)

        % Simply assign the rate
        r1_rates(n)=current_r1_rate;

    % 3x3 tensor specification
    elseif isnumeric(current_r1_rate)&&...
           (size(current_r1_rate,1)==3)&&...
           (size(current_r1_rate,2)==3)

        % Make sure the angles are specified
        if ~exist('euler_angles','var')
            error('Euler angles must be specified with anisotropic T1/T2 relaxation theory.');
        end

        % Compute orientation ort (this matches alphas=0 of two-angle grids)
        ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3));

        % Get the rate at the current orientation
        r1_rates(n)=ort*current_r1_rate*ort';

    % Function handle specification
    elseif isa(current_r1_rate,'function_handle')

        % Make sure the angles are specified
        if ~exist('euler_angles','var')
            error('Euler angles must be specified with anisotropic T1/T2 relaxation theory.');
        end

        % Call the function handle
        r1_rates(n)=current_r1_rate(euler_angles(1),euler_angles(2),euler_angles(3));

    else

        % Complain and bomb out
        error('unknown R1 rate specification.');

    end

end

% Fill in R2 relaxation rates
for n=1:numel(spin_system.comp.isotopes)

    % Get the specification
    current_r2_rate=spin_system.rlx.r2_rates{n};

    % Scalar relaxation rate
    if isnumeric(current_r2_rate)&&isscalar(current_r2_rate)

        % Simply assign the rate
        r2_rates(n)=current_r2_rate;

    % 3x3 tensor specification
    elseif isnumeric(current_r2_rate)&&...
           (size(current_r2_rate,1)==3)&&...
           (size(current_r2_rate,2)==3)

        % Make sure the angles are specified
        if ~exist('euler_angles','var')
            error('Euler angles must be specified with anisotropic T1/T2 relaxation theory.');
        end

        % Compute orientation ort (this matches alphas=0 of two-angle grids)
        ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3));

        % Get the rate at the current orientation
        r2_rates(n)=ort*current_r2_rate*ort';

    % Function handle specification
    elseif isa(current_r2_rate,'function_handle')

        % Make sure the angles are specified
        if ~exist('euler_angles','var')
            error('Euler angles must be specified with anisotropic T1/T2 relaxation theory.');
        end

        % Call the function handle
        r2_rates(n)=current_r2_rate(euler_angles(1),euler_angles(2),euler_angles(3));

    else

        % Complain and bomb out
        error('unknown R1 rate specification.');

    end

end

% Make sure the rates make sense
if any(~isreal(r1_rates),'all')||any(~isreal(r2_rates),'all')
    error('all R1 and R2 relaxation rates must be real numbers.');
end

% Preallocate direct-sum relaxation diagonals
matrix_dim=spin_system.bas.offsets(end);
r1_diagonal=zeros(matrix_dim,1);
r2_diagonal=zeros(matrix_dim,1);

% Sum the local single-spin rates within each substance
for n=1:spin_system.bas.nsubst
    [L,M]=lin2lm(spin_system.bas.basis{n});
    spins=spin_system.chem.parts{n};
    idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    r1_diagonal(idx)=((L~=0)&(M==0))*r1_rates(spins,1);
    r2_diagonal(idx)=((L~=0)&(M~=0))*r2_rates(spins,1);
end

% Build relaxation superoperators
R1Op=-spdiags(r1_diagonal,0,matrix_dim,matrix_dim);
R2Op=-spdiags(r2_diagonal,0,matrix_dim,matrix_dim);

end

% Consistency enforcement
function grumble(spin_system)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available in sphten-liouv formalism.');
end
end

% The state that separates its scholars from its warriors will
% have its thinking done by cowards and its fighting by fools.
%
% Sir William Francis Bacon

