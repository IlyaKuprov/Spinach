% Steady-state saturated continuous-wave detection - calculates the
% expectation value of the detection state in the steady state of a
% thermalised spin system driven by a continuous irradiation field
% that is retained to all orders in its amplitude. Unlike slowpass.m,
% which returns a linear response resolvent matrix element and is
% therefore first order in the irradiation amplitude by construction,
% this function solves, once per point of the offset grid, the inhomo-
% geneous equation obeyed by the steady-state density operator under
% the drive. Syntax:
%
%         spectrum=slowpass_sat(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep         vector with two elements giving
%                             the spectrum frequency extents
%                             in Hz
%
%    parameters.npoints       number of points in the spectrum
%
%    parameters.spins         cell array with one character
%                             string, specifying the isotope
%                             that the irradiation acts upon
%
%    parameters.irr_powers    irradiation nutation frequency,
%                             in Hz
%
%    parameters.coil          detection state
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    spectrum  - a complex column vector of length parameters.npoints,
%                containing the steady-state expectation value of the
%                detection state at each point of the offset frequency
%                interval requested in parameters.sweep
%
% Note: the signature matches slowpass.m, but this solver requires
%       spherical-tensor Liouville space, thermalised relaxation, one
%       irradiated isotope, and a single detection column. The scalar
%       irr_powers specifies nu_1; the drive is 2*pi*nu_1*Lx. rho0 is
%       ignored: the source is supplied by thermalised relaxation.
%
% Note: the amplitude convention differs from slowpass.m, which scales
%       its answer to match the unnormalised FFT of a free induction
%       decay. A steady-state expectation value needs no such scaling
%       and none is applied here.
%
% Note: the relaxation superoperator must be thermalised, meaning that
%       inter.equilibrium must be set to 'IME' or 'dibari'. With the
%       Spinach default of 'zero' the steady state of a driven dissi-
%       pative system is identically zero; this function detects that
%       situation and refuses to run.
%
% Note: retain the longitudinal and both transverse coherence orders
%       in the basis. This function performs no trajectory reduction;
%       a basis truncated for linear response need not support saturation.
%
% Note: spins without a relaxation pathway produce exactly zero pivots
%       in the longitudinal block, which stops ilu; the preconditioner
%       is therefore built from a slightly shifted matrix. The shift
%       does not enter the matrix that is solved.
%
% Note: thermal equilibrium is assumed to be isotropic, meaning that
%       it occupies the rank zero coordinate of the spatial basis. An
%       orienting potential would place equilibrium population in the
%       higher ranks as well, and is not supported here. Spatial inputs
%       must use the gridfree space-times-spin ordering, with the uniform
%       distribution first. Other spatial contexts are not supported.
%
% Note: every substance concentration must be individually conserved;
%       population-changing kinetics are refused. After these conserved
%       directions are pinned, the steady state must be unique.
%
% Note: explicit matrices are required. CPU backslash is used up to
%       100000 active coordinates; larger problems use shifted-ILU
%       right-preconditioned GMRES with tolerance 1e-6 and restart length
%       10. ILU factors are recomputed at each frequency because the sweep
%       generator is not the identity and leaves populations unshifted.
%       The equation is scaled by its median nonzero diagonal magnitude;
%       the unpreconditioned relative residual must also be at most 1e-6.
%       GPU and polyadic solver paths are not implemented.
%
% Reference: https://doi.org/10.1007/s00723-021-01390-7
%
% alexey.bogdanov@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function spectrum=slowpass_sat(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Get the offset frequency grid in angular units
freq_grid=2*pi*linspace(parameters.sweep(1),...
                        parameters.sweep(2),...
                        parameters.npoints)';

% Build the irradiation and the sweep superoperators in spin space
irr_op=operator(spin_system,'Lx',parameters.spins{1});
sweep_op=operator(spin_system,'Lz',parameters.spins{1});

% Embed the spin superoperators in every spatial basis coordinate
if isfield(parameters,'spc_dim')
    irr_op=kron(speye(parameters.spc_dim),irr_op);
    sweep_op=kron(speye(parameters.spc_dim),sweep_op);
end

% Compose the part of the Liouvillian that does not depend on the sweep
L=H+1i*R+1i*K+2*pi*parameters.irr_powers*irr_op;

% Locate the rank-zero unit coordinate of every substance
unit_idx=spin_system.bas.offsets(1:end-1)+1;

% Require a thermal source
if norm(L(:,unit_idx),1)<=spin_system.tols.liouv_zero
    error('relaxation superoperator is not thermalised, set inter.equilibrium.');
end

% Refuse to run when the unit coordinate is not a conserved quantity
if norm(L(unit_idx,:),1)>spin_system.tols.liouv_zero
    error('the unit coordinate is not conserved, check relaxation and kinetics.');
end

% Pin the conserved substance concentrations
active=true(size(L,1),1); active(unit_idx)=false;

% Thermal source term, independent of the sweep and of the irradiation
source=-L(active,unit_idx)*spin_system.chem.concs(:);

% Project the Liouvillian, the sweep generator, and the detection state
L=L(active,active); sweep_op=sweep_op(active,active);
coil=parameters.coil(active);

% Include the known contribution of the pinned coordinates
spectrum=(parameters.coil(unit_idx)'*spin_system.chem.concs(:))*ones(size(freq_grid));

% Use the iterative solver for large active spaces
if size(L,1)>1e5

    % Inform the user
    report(spin_system,'using preconditioned CPU GMRES path...');

    % Get the scale of the problem from the nonzero diagonal elements
    diag_mags=full(abs(diag(L))); diag_scale=median(diag_mags(diag_mags>0));

    % Scale the equation before iterative solution
    L=L/diag_scale; sweep_op=sweep_op/diag_scale; source=source/diag_scale;

    % Factor each shifted field point and measure the physical residual
    parfor n=1:numel(freq_grid)
        A=L+freq_grid(n)*sweep_op;
        [M1,M2]=ilu(A+1e-6*speye(size(A)),struct('type','crout','droptol',1e-3));
        action=@(x)A*(M2\(M1\x));
        [rho,flag]=gmres(action,full(source),10,1e-6,numel(source));
        rho=M2\(M1\rho);
        residual=norm(A*rho-source);
        if (flag~=0)||(residual>1e-6*norm(source))
            error('steady-state GMRES did not converge at offset %g Hz.',freq_grid(n)/(2*pi));
        end
        spectrum(n)=spectrum(n)+coil'*rho;
    end

else

    % Inform the user
    report(spin_system,'using CPU backslash path...');

    % Run the sweep in parallel
    parfor n=1:numel(freq_grid)
        spectrum(n)=spectrum(n)+coil'*((L+freq_grid(n)*sweep_op)\source);
    end

end

% Refuse a non-finite stationary solution
if any(~isfinite(spectrum)), error('steady-state solve returned non-finite values.'); end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (size(H,1)~=size(H,2))||(~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must be square and of the same dimension.');
end
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function requires bas.formalism=''sphten-liouv''.');
end
if ~isfield(parameters,'sweep')
    error('spectral range must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||(numel(parameters.sweep)~=2)||...
   any(~isfinite(parameters.sweep))
    error('parameters.sweep vector must have two finite real elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<2)||(mod(parameters.npoints,1)~=0)||...
   (~isfinite(parameters.npoints))
    error('parameters.npoints should be an integer greater than one.');
end
if ~isfield(parameters,'spins')
    error('the irradiated isotope must be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||(~ischar(parameters.spins{1}))
    error('parameters.spins must be a cell array with one character string.');
end
if ~isfield(parameters,'irr_powers')
    error('the nutation frequency must be specified in parameters.irr_powers variable.');
end
if (~isnumeric(parameters.irr_powers))||(~isreal(parameters.irr_powers))||...
   (~isscalar(parameters.irr_powers))||(parameters.irr_powers<0)||(~isfinite(parameters.irr_powers))
    error('parameters.irr_powers must be a finite non-negative real scalar.');
end
if ~isfield(parameters,'coil')
    error('the detection state must be specified in parameters.coil variable.');
end
if (~isnumeric(parameters.coil))||(~iscolumn(parameters.coil))||(size(parameters.coil,1)~=size(H,1))||...
   any(~isfinite(parameters.coil))
    error('parameters.coil must be a finite numeric column with as many rows as H.');
end
end


