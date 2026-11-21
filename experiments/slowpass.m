% Slow passage detection - calculates spectrum values at the user-
% specified frequency positions using the Fourier transform of the
% Liouville - von Neumann equation. The biggest advantage over the
% fid+fft style detection is easy parallelization and the possibi-
% lity of getting spectrum values at specific frequencies without
% recalculating the entire free induction decay. Syntax:
%
%         spectrum=slowpass(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep         vector with two elements giving
%                             the spectrum frequency extents
%                             in Hz
%
%    parameters.npoints       number of points in the spectrum
%
%    parameters.mode          'linear' (default) or 'saturated'
%
%    parameters.spins         saturated mode: one irradiated isotope
%                             in a cell array
%
%    parameters.irr_powers    saturated mode: nutation frequency in Hz
%
%    parameters.rho0          initial state (linear mode only)
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
%    spectrum  - the spectrum of the system with the specified
%                starting state detected on the specified coil
%                state within the frequency interval requested
%
% Saturated mode returns the stationary detector expectation without FFT
% scaling. Use thermalised relaxation, explicit spherical-tensor Liouville
% matrices, and a basis retaining populations and both transverse orders.
% Substance populations must be conserved. For gridfree the equilibrium
% distribution must be uniform. See doi:10.1007/s00723-021-01390-7.
%
% Linear mode: relaxation must be present in the system dynamics, or the
%       matrix inversion operation would fail to converge. The re-
%       laxation matrix R must *not* be thermalized.
%
% Note: Liouville-space identity components are excluded only when
%       the identity sector is decoupled from spin order in both
%       directions, to within tols.liouv_zero. Its resolvent block
%       is shifted by 1 inverse second to remove stationary poles;
%       the spin-order block is unchanged. Spatial contexts supply
%       parameters.spc_dim for the space-times-spin embedding.
%       Coupled identity sectors, including selective reactions,
%       and wavefunction inputs retain their original resolvent.
%
% alexey.bogdanov@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=slowpass.m>

function spectrum=slowpass(spin_system,parameters,H,R,K)

% Check consistency
if ~isfield(parameters,'mode'), parameters.mode='linear'; end
grumble('inputs');

% Get the frequency grid
freq_grid=2*pi*linspace(parameters.sweep(1),...
                        parameters.sweep(2),...
                        parameters.npoints)';

% Prepare the offset-independent quantities for each experiment
saturated=strcmp(parameters.mode,'saturated');
if saturated

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

    % Pin the conserved substance concentrations
    active=true(size(L,1),1);
    active(unit_idx)=false;

    % Thermal source term, independent of the sweep and of the irradiation
    rho0_subs=-L(active,unit_idx)*spin_system.chem.concs(:);

    % Project the Liouvillian, the sweep generator, and the detection state
    L_subs=L(active,active);
    sweep_op=sweep_op(active,active);
    coil_subs=parameters.coil(active);

    % Include the known contribution of the pinned coordinates
    spectrum=(parameters.coil(unit_idx)'*spin_system.chem.concs(:))*ones(size(freq_grid));

    nblocks=1;

else

    % Move into adjoint representation if needed
    [spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K);

    % Compose the Liouvillian
    L=H+1i*R+1i*K;

    % Compute subspace projectors
    projectors=reduce(spin_system,L,parameters.coil);

    % Get identity directions only in Liouville-space formalisms
    U=sparse(size(L,1),0);
    if ~strcmp(spin_system.bas.formalism,'zeeman-wavef')

        % Get normalised unit states, one per substance
        U=cell(1,spin_system.bas.nsubst);
        for n=1:spin_system.bas.nsubst
            unit_system=spin_system; unit_system.chem.concs(:)=0;
            unit_system.chem.concs(n)=1; U{n}=unit_state(unit_system);
        end
        U=[U{:}];

        % Embed the spin identities in every spatial basis coordinate
        if isfield(parameters,'spc_dim')
            U=kron(speye(parameters.spc_dim),U);
        end
    end

    spectrum=zeros(size(freq_grid),'like',1i);
    nblocks=numel(projectors);

end

% Solve the same shifted system in every contributing block
for k=1:nblocks

    % Linear response uses a projected resolvent in each detection subspace
    if ~saturated

        % Project into the current subspace
        rho0_subs=projectors{k}'*parameters.rho0;
        coil_subs=projectors{k}'*parameters.coil;
        L_subs=projectors{k}'*L*projectors{k};

        % Normalise the projected identity directions
        U_subs=projectors{k}'*U; U_subs=U_subs(:,any(U_subs,1));
        U_subs=U_subs./sqrt(sum(abs(U_subs).^2,1));

        % Remove identities only when they decouple from spin order in both directions
        unit_block=U_subs'*L_subs*U_subs;
        if norm(L_subs*U_subs-U_subs*unit_block,1)<=spin_system.tols.liouv_zero&&...
           norm(U_subs'*L_subs-unit_block*U_subs',1)<=spin_system.tols.liouv_zero
            rho0_subs=rho0_subs-U_subs*(U_subs'*rho0_subs);
            coil_subs=coil_subs-U_subs*(U_subs'*coil_subs);
            L_subs=L_subs-1i*(U_subs*U_subs');
        end

        L_subs=1i*L_subs;
        sweep_op=1i*speye(size(L_subs));

    end

    % Select the backend before sweeping the offsets
    use_gpu=~saturated&&ismember('gpu',spin_system.sys.enable);
    iterative=~use_gpu&&(size(rho0_subs,2)==1)&&...
              ((saturated&&(size(L_subs,1)>1e5))||(~saturated&&(size(L_subs,1)>5000)));
    M1=[]; M2=[];
    if use_gpu
        report(spin_system,'using GPU backslash path...');
        L_subs=gpuArray(full(L_subs)); sweep_op=gpuArray(sweep_op);
        rho0_subs=gpuArray(full(rho0_subs)); coil_subs=gpuArray(full(coil_subs));
    elseif iterative
        report(spin_system,'using preconditioned CPU GMRES path...');
        if saturated
            diag_mags=full(abs(diag(L_subs)));
            diag_scale=median(diag_mags(diag_mags>0));
            L_subs=L_subs/diag_scale; sweep_op=sweep_op/diag_scale;
            rho0_subs=full(rho0_subs/diag_scale);
        else
            [M1,M2]=ilu(L_subs+mean(freq_grid)*sweep_op,struct('type','crout','droptol',1e-3));
        end
    else
        report(spin_system,'using CPU backslash path...');
    end

    % Solve, detect, and accumulate each frequency point
    flags=zeros(size(freq_grid)); residuals=zeros(size(freq_grid));
    parfor n=1:numel(freq_grid)
        A=L_subs+freq_grid(n)*sweep_op;
        if iterative
            if saturated
                [P,Q]=ilu(A+1e-6*speye(size(A)),struct('type','crout','droptol',1e-3));
                action=@(x)A*(Q\(P\x));
            else
                action=A;
            end
            [rho,flag]=gmres(action,rho0_subs,10,1e-6,numel(rho0_subs),M1,M2);
            flags(n)=flag;
            if saturated
                rho=Q\(P\rho);
                residuals(n)=norm(A*rho-rho0_subs);
            end
        else
            rho=A\rho0_subs;
        end
        response=dot(coil_subs,rho);
        if use_gpu, response=gather(response); end
        spectrum(n)=spectrum(n)+response;
    end
    grumble('solution');

end

% Only the linear resolvent uses the unnormalised FFT convention
if ~saturated
    spectrum=spectrum*(abs(parameters.sweep(2)-parameters.sweep(1))*...
                       parameters.npoints/(parameters.npoints-1));
end

% Consistency enforcement
    function grumble(stage)
        if strcmp(stage,'solution')
            if ~saturated
                failed=find(flags~=0,1);
                if ~isempty(failed)
                    warning('linear GMRES returned flag %d at offset %g Hz.',flags(failed),freq_grid(failed)/(2*pi));
                end
                return;
            end
            failed=find((flags~=0)|(residuals>1e-6*norm(rho0_subs)),1);
            if ~isempty(failed)
                error('steady-state GMRES did not converge at offset %g Hz.',freq_grid(failed)/(2*pi));
            end
            if any(~isfinite(spectrum)), error('steady-state solve returned non-finite values.'); end
            return;
        end
        if (~ischar(parameters.mode))||(~isrow(parameters.mode))||...
           (~ismember(parameters.mode,{'linear','saturated'}))
            error('parameters.mode must be either linear or saturated.');
        end
        if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
           (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
            error('H, R and K arguments must be matrices.');
        end
        if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
            error('H, R and K matrices must have the same dimension.');
        end
        if ~isfield(parameters,'sweep')
            error('spectral range must be specified in parameters.sweep variable.');
        end
        if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||(numel(parameters.sweep)~=2)
            error('parameters.sweep vector must have two real elements');
        end
        if ~isfield(parameters,'npoints')
            error('number of points must be specified in parameters.npoints variable.');
        end
        if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
           (~isreal(parameters.npoints))||(parameters.npoints<2)||(mod(parameters.npoints,1)~=0)
            error('parameters.npoints should be an integer greater than one.');
        end
        if ~isfield(parameters,'coil')
            error('the detection state must be specified in parameters.coil variable.');
        end
        if (~isnumeric(parameters.coil))||(~ismatrix(parameters.coil))
            error('parameters.coil must be a numeric matrix.');
        end
        if size(parameters.coil,1)~=size(H,2)
            error('parameters.coil must have the same number of rows as H.');
        end
        if strcmp(parameters.mode,'linear')
            if ~isfield(parameters,'rho0')
                error('the initial state must be specified in parameters.rho0 variable.');
            end
            if (~isnumeric(parameters.rho0))||(~ismatrix(parameters.rho0))
                error('parameters.rho0 must be a numeric matrix.');
            end
            if size(parameters.rho0,1)~=size(H,2)
                error('parameters.rho0 must have the same number of rows as H.');
            end
        else
            if size(H,1)~=size(H,2)
                error('saturated mode requires square generators.');
            end
            if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
                error('saturated mode requires spherical-tensor Liouville space.');
            end
            if any(~isfinite(parameters.sweep))
                error('parameters.sweep must contain finite offsets.');
            end
            if (~iscolumn(parameters.coil))||any(~isfinite(nonzeros(parameters.coil)))
                error('saturated mode requires a single finite detection column.');
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
            unit_idx=spin_system.bas.offsets(1:end-1)+1;
            if norm(H(:,unit_idx)+1i*R(:,unit_idx)+1i*K(:,unit_idx),1)<=spin_system.tols.liouv_zero
                error('relaxation superoperator is not thermalised, set inter.equilibrium.');
            end
            if norm(H(unit_idx,:)+1i*R(unit_idx,:)+1i*K(unit_idx,:),1)>spin_system.tols.liouv_zero
                error('the unit coordinate is not conserved, check relaxation and kinetics.');
            end
        end
    end

end

% "Our first rule here, Miss Taggart," he answered, "is that
%  one must always see for oneself."
%
% Ayn Rand, "Atlas Shrugged"


