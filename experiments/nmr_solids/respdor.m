% R- and S-RESPDOR echoes for the Liouville-space singlerot context.
% The R variant uses a two-rotor-period rotary-resonance block;
% the S variant uses a six-period phase-supercycled R4 block.
% The two recoupling arms have the same pulse order. A requested
% zero count returns the prepared-state overlap without an echo.
%
% Syntax:
%
%             curve=respdor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins     - two distinct isotope labels in a cell row,
%                           observed spin-1/2 channel first
%    parameters.sequence  - 'R' or 'S'
%    parameters.rate      - non-zero signed MAS rate in Hz
%    parameters.ncycles   - non-empty row of non-negative integer
%                           recoupling block counts, per arm
%    parameters.echo_rf   - positive observed-channel pi-pulse RF in Hz
%    parameters.sat_rf    - non-negative quadrupolar irradiation RF in Hz
%    parameters.sat_time  - irradiation window in seconds, at least
%                           the pi-pulse length and at most two periods
%    parameters.spc_dim   - spatial dimension supplied by singlerot
%    parameters.rho0      - initial column state, lifted by singlerot
%    parameters.coil      - detection column state, lifted by singlerot
%    H,R,K                - explicit lifted drift, relaxation, and
%                           kinetics matrices supplied by singlerot
%
% Outputs:
%
%    curve                - complex rows: reference, irradiated echo,
%                           and reference minus irradiated echo
%
% talos@spindynamics.org

function curve=respdor(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Lift the RF generators into rotor space
L=H+1i*R+1i*K;
space_id=speye(parameters.spc_dim);
Ix=kron(space_id,operator(spin_system,'Lx',parameters.spins{1}));
Iy=kron(space_id,operator(spin_system,'Ly',parameters.spins{1}));
Sx=kron(space_id,operator(spin_system,'Lx',parameters.spins{2}));

% Set the rotor-synchronised recoupling block
period=1/abs(parameters.rate);
if strcmp(parameters.sequence,'R')

    % Two periods of rotary resonance along X
    phases=0;
    block_rf=1/period;
    dwell=2*period;

else

    % Three phase increments of the R4 block and its phase inversion
    phases=deg2rad([90 -90 90 -90 -90 90 -90 90].'+[0 120 240]);
    block_rf=2/period;
    dwell=period/4;

end
block_gen=cell(1,numel(phases));
for n=1:numel(phases)

    % Keep the full drift during every recoupling pulse
    block_gen{n}=L+2*pi*block_rf*(cos(phases(n))*Ix+sin(phases(n))*Iy);

end

% Centre the observed-channel pi pulse in the common irradiation window
echo_width=0.5/parameters.echo_rf;
outer_delay=period-parameters.sat_time/2;
sat_flank=(parameters.sat_time-echo_width)/2;
L_sat=L+2*pi*parameters.sat_rf*Sx;
L_echo=L+2*pi*parameters.echo_rf*Ix;
L_both=L_sat+2*pi*parameters.echo_rf*Ix;

% Advance the state and back-propagate the receiver by one block per count
rho=parameters.rho0;
coil=parameters.coil;
echoes=repmat(coil'*rho,2,max(parameters.ncycles)+1);
for n=1:max(parameters.ncycles)

    % Apply the same chronological recoupling block in both arms
    for k=1:numel(block_gen)
        rho=step(spin_system,block_gen{k},rho,dwell);
        coil=step(spin_system,block_gen{end+1-k}',coil,-dwell);
    end

    % Refocus both branches with identical elapsed time
    pre=step(spin_system,L,rho,outer_delay);
    ref=step(spin_system,L,pre,sat_flank);
    ref=step(spin_system,L_echo,ref,echo_width);
    ref=step(spin_system,L,ref,sat_flank);
    deph=step(spin_system,L_sat,pre,sat_flank);
    deph=step(spin_system,L_both,deph,echo_width);
    deph=step(spin_system,L_sat,deph,sat_flank);
    post=step(spin_system,L,[ref deph],outer_delay);
    echoes(:,n+1)=(coil'*post).';

end

% Select requested counts and retain the signed complex difference
echoes=echoes(:,parameters.ncycles+1);
curve=[echoes; echoes(1,:)-echoes(2,:)];

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('RESPDOR requires sphten-liouv or zeeman-liouv formalism.');
end
if ~isnumeric(parameters.rho0)||~iscolumn(parameters.rho0)||...
   ~isnumeric(parameters.coil)||~isequal(size(parameters.coil),size(parameters.rho0))
    error('parameters.rho0 and parameters.coil must be matching column vectors.');
end
if ~isnumeric(H)||~ismatrix(H)||~isnumeric(R)||~isnumeric(K)||...
   ~isequal(size(H),size(R),size(K),[numel(parameters.rho0) numel(parameters.rho0)])
    error('H, R, and K must be explicit matrices matching the lifted state.');
end
if ~iscell(parameters.spins)||~isequal(size(parameters.spins),[1 2])||...
   ~all(cellfun(@ischar,parameters.spins))||...
   any(~ismember(parameters.spins,spin_system.comp.isotopes))||...
   strcmp(parameters.spins{1},parameters.spins{2})
    error('parameters.spins must contain two distinct isotope labels present in the system.');
end
if any(spin_system.comp.mults(strcmp(spin_system.comp.isotopes,parameters.spins{1}))~=2)||...
   any(spin_system.comp.mults(strcmp(spin_system.comp.isotopes,parameters.spins{2}))<3)
    error('RESPDOR requires spin 1/2 on the first channel and a quadrupolar spin on the second.');
end
if ~ischar(parameters.sequence)||~ismember(parameters.sequence,{'R','S'})
    error('parameters.sequence must be ''R'' or ''S''.');
end
validateattributes(parameters.rate,{'numeric'},{'scalar','real','finite','nonzero'});
validateattributes(parameters.ncycles,{'numeric'},{'row','nonempty','real','finite','nonnegative','integer'});
validateattributes(parameters.echo_rf,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(parameters.sat_rf,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(parameters.sat_time,{'numeric'},{'scalar','real','finite'});
if parameters.sat_time<0.5/parameters.echo_rf||parameters.sat_time>2/abs(parameters.rate)
    error('parameters.sat_time must be between the pi-pulse length and two rotor periods.');
end
end


