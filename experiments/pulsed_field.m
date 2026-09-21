% Magnetisation dynamics under a time-dependent magnetic field along
% the Z axis of the laboratory frame with spin-phonon relaxation, as
% measured in pulsed-field magnetometry of molecular magnets. The
% field profile is replaced by a staircase; on each stair the Hamil-
% tonian is constant, the spin-phonon dissipator is rebuilt in the
% eigenbasis of that Hamiltonian, and the density matrix is propa-
% gated in that eigenbasis by a symmetric split: exact coherent
% phases for half a stair, the dissipative step to second order in
% the dissipator times the stair width, and the phases again. The
% dissipator times the stair width must be small; the coherent part
% is treated exactly for any stair width. The dissipator is applied
% as Hilbert space matrix products (see phonon_oper.m), so the cost
% of a stair is cubic in the dimension of the Hilbert space. Syntax:
%
%          answer=pulsed_field(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.field_prof - function handle returning the field
%                            in Tesla at a time in seconds
%
%    parameters.hzeeman    - Zeeman operator per Tesla, rad/s/T,
%                            Hilbert space, supplied by the con-
%                            text when 'zeeman_op' is requested
%                            in parameters.needs
%
%    parameters.timestep   - stair width, seconds
%
%    parameters.nsteps     - number of stairs
%
%    parameters.coil       - Hermitian Hilbert space observable
%                            operator or a cell array of them
%
%    parameters.phonon_x   - spin-phonon coupling operator, see
%                            rlx_phonon.m
%
%    parameters.phonon_i0  - phonon spectral density prefactor,
%                            see rlx_phonon.m
%
%    parameters.phonon_alpha - phonon spectral density exponent,
%                              1 or above, see rlx_phonon.m
%
%    parameters.nout       - number of stairs between recorded
%                            observable values
%
%    H - Hamiltonian received from the context function, Hilbert
%        space, containing the Zeeman term at sys.magnet=1 Tesla;
%        the function removes that term and adds the field of
%        each stair itself
%
%    R - relaxation superoperator received from the context
%        function; ignored, the spin-phonon dissipator is built
%        here at every stair
%
%    K - kinetics superoperator received from the context
%        function; ignored
%
% Outputs:
%
%    answer.t     - column of recording times, seconds
%
%    answer.field - column of field values at those times, Tesla
%
%    answer.obs   - matrix of observable expectation values, one
%                   column per coil, at the recording times
%
% Note: the sequence works in zeeman-hilb formalism under the crystal
%       and powder contexts, which assemble the anisotropic part of
%       the Hamiltonian; the liquid context drops that part, and with
%       it the crystal field of a giant spin. The context must be
%       called with the labframe assumption set, so that H and the
%       Zeeman operator are built consistently; the powder context
%       must be called with parameters.sum_up=false because the
%       answer is a structure; additional rotating frames (parame-
%       ters.rframes) and frequency offsets (parameters.offset) are
%       not supported because the field operator is added in the
%       laboratory frame. The temperature of the phonon bath is
%       inter.temperature.
%
% Note: sys.magnet must be 1 Tesla, so that parameters.hzeeman is
%       the Zeeman operator per Tesla; the Hamiltonian received from
%       the context then contains the Zeeman term at 1 Tesla, which
%       this function removes before adding the field on each stair.
%       The initial state is the thermal equilibrium of the field-
%       free Hamiltonian at the temperature of the phonon bath.
%
% Note: the Hamiltonian on each stair uses the field at the midpoint
%       of the stair; answer.field is the profile evaluated at the
%       recording times, which are the ends of the recorded stairs.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pulsed_field.m>

function answer=pulsed_field(spin_system,parameters,H,R,K) %#ok<INUSD>

% Check consistency
grumble(spin_system,parameters,H);

% Put the coils into a cell array
if iscell(parameters.coil), coils=parameters.coil; else, coils={parameters.coil}; end

% Preallocate the output
nrec=floor(parameters.nsteps/parameters.nout);
answer.t=zeros(nrec,1); answer.field=zeros(nrec,1); answer.obs=zeros(nrec,numel(coils));

% Remove the unit field Zeeman term supplied by the context
H=H-parameters.hzeeman; H=(H+H')/2;

% Thermal equilibrium at zero field as the initial state
[V,E]=eig(full(H),'vector'); pops=exp(-spin_system.tols.hbar*(E-min(E))/(spin_system.tols.kbol*spin_system.rlx.temperature));
rho=V*diag(pops/sum(pops))*V'; dt=parameters.timestep; nrec=0;

% Loop over the stairs
for n=1:parameters.nsteps

    % Field at the stair midpoint and the Hamiltonian on the stair
    field=parameters.field_prof((n-0.5)*dt);
    if (~isnumeric(field))||(~isreal(field))||(~isscalar(field))||(~isfinite(field))
        error('parameters.field_prof must return a real finite scalar.');
    end
    H_curr=H+field*parameters.hzeeman; H_curr=full((H_curr+H_curr')/2);

    % Eigensystem of the stair Hamiltonian and the coherent half-stair phases
    [V,E]=eig(H_curr,'vector'); phases=exp(-1i*(E-E.')*dt/2);

    % Spin-phonon coupling operator and its thermally dressed form in the eigenbasis
    XE=V'*parameters.phonon_x*V; XE=(XE+XE')/2;
    RE=phonon_oper(spin_system,E,XE,parameters.phonon_i0,parameters.phonon_alpha,spin_system.rlx.temperature);

    % Dissipator as matrix products in the eigenbasis
    dissip=@(rho)-pi*(XE*(RE*rho)-(RE*rho)*XE+rho*(RE'*XE)-(XE*rho)*RE');

    % Symmetric split step in the eigenbasis
    rho_eig=phases.*(V'*rho*V); drho=dissip(rho_eig);
    rho_eig=rho_eig+dt*drho+(dt^2/2)*dissip(drho);
    rho=V*(phases.*rho_eig)*V';

    % Record the field at the end of the stair, the observables, and the progress
    if mod(n,parameters.nout)==0
        nrec=nrec+1; answer.t(nrec)=n*dt; field=parameters.field_prof(n*dt);
        if (~isnumeric(field))||(~isreal(field))||(~isscalar(field))||(~isfinite(field))
            error('parameters.field_prof must return a real finite scalar.');
        end
        answer.field(nrec)=field;
        for k=1:numel(coils)
            answer.obs(nrec,k)=real(trace(coils{k}'*rho));
        end
        report(spin_system,['stair ' num2str(n) ' of ' num2str(parameters.nsteps) ', field ' ...
                            num2str(field) ' T, observable ' num2str(answer.obs(nrec,1))]);
    end

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('this function is only available in zeeman-hilb formalism.');
end
if spin_system.inter.magnet~=1
    error('sys.magnet must be 1 Tesla, the field is set by parameters.field_prof.');
end
if ~strcmp(spin_system.inter.assumptions,'labframe')
    error('this function requires the labframe assumption set in the context call.');
end
if (~isnumeric(H))||(size(H,1)~=size(H,2))
    error('H must be a square matrix.');
end
if ~isfield(parameters,'field_prof')||(~isa(parameters.field_prof,'function_handle'))
    error('parameters.field_prof must be a function handle returning the field in Tesla.');
end
if ~isfield(parameters,'hzeeman')||(~isnumeric(parameters.hzeeman))||any(size(parameters.hzeeman)~=size(H))
    error('parameters.hzeeman must be a matrix of the same dimension as H, add ''zeeman_op'' to parameters.needs.');
end
if ~isfield(parameters,'timestep')||(~isnumeric(parameters.timestep))||(~isreal(parameters.timestep))||...
   (~isscalar(parameters.timestep))||(~isfinite(parameters.timestep))||(parameters.timestep<=0)
    error('parameters.timestep must be a positive real scalar.');
end
if ~isfield(parameters,'nsteps')||(~isnumeric(parameters.nsteps))||(~isscalar(parameters.nsteps))||(mod(parameters.nsteps,1)~=0)||(parameters.nsteps<1)
    error('parameters.nsteps must be a positive integer.');
end
if ~isfield(parameters,'nout')||(~isnumeric(parameters.nout))||(~isscalar(parameters.nout))||(mod(parameters.nout,1)~=0)||(parameters.nout<1)
    error('parameters.nout must be a positive integer.');
end
if ~isfield(parameters,'coil')||(~(isnumeric(parameters.coil)||iscell(parameters.coil)))
    error('parameters.coil must be an observable operator or a cell array of them.');
end
if iscell(parameters.coil), coils=parameters.coil; else, coils={parameters.coil}; end
if isempty(coils)
    error('parameters.coil must contain at least one observable operator.');
end
for k=1:numel(coils)
    if (~isnumeric(coils{k}))||any(size(coils{k})~=size(H))||any(~isfinite(coils{k}(:)))||(~ishermitian(coils{k}))
        error('every coil must be a Hermitian matrix of the same dimension as H with finite elements.');
    end
end
if isfield(parameters,'rframes')&&(~isempty(parameters.rframes))
    error('additional rotating frames are not supported by this function.');
end
if isfield(parameters,'offset')&&any(parameters.offset(:)~=0)
    error('frequency offsets are not supported by this function.');
end
if isfield(parameters,'sum_up')&&parameters.sum_up
    error('the powder context must be called with parameters.sum_up=false.');
end
if ~isfield(parameters,'phonon_x')||(~isnumeric(parameters.phonon_x))||any(size(parameters.phonon_x)~=size(H))||(~ishermitian(parameters.phonon_x))
    error('parameters.phonon_x must be a Hermitian matrix of the same dimension as H.');
end
if ~isfield(parameters,'phonon_i0')||(~isnumeric(parameters.phonon_i0))||(~isscalar(parameters.phonon_i0))||(parameters.phonon_i0<0)
    error('parameters.phonon_i0 must be a non-negative real scalar.');
end
if ~isfield(parameters,'phonon_alpha')||(~isnumeric(parameters.phonon_alpha))||(~isscalar(parameters.phonon_alpha))||(parameters.phonon_alpha<1)
    error('parameters.phonon_alpha must be a real scalar not smaller than 1.');
end
if ~isfield(spin_system.rlx,'temperature')||isempty(spin_system.rlx.temperature)
    error('the phonon bath temperature must be specified in inter.temperature.');
end
end

% Experience is the name everyone gives to their mistakes.
%
% Oscar Wilde

