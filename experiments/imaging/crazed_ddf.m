% Ideal two-pulse CRAZED experiment with a self-consistent secular
% distant dipolar field in a periodic one-dimensional sample. Syntax:
%
%              fid=imaging(spin_system,@crazed_ddf,parameters)
%
% Parameters (in addition to the ordinary imaging parameters):
%
%   spins          - one isotope label in a cell; all spins must match
%   direction      - unit three-vector along the grating in lab axes,
%                    with the static magnetic field along Z
%   demag_factor   - axial demagnetising factor in [0,1]
%   magnetisation  - equilibrium magnetisation density in A/m;
%                    zero switches off feedback without changing rho0
%   rho_ref        - local equilibrium state column corresponding to
%                    that density, with a real, nonzero Z projection
%   flip_angles    - two signed X rotation angles in radians
%   gradient_areas - two signed instantaneous gradient areas in T s/m;
%                    each must wind an integer number of pitches
%   sweep          - acquisition sweep width in Hz
%   npoints        - number of complex samples, including time zero
%   ddf_steps      - positive integer, LG4 substeps per sampling interval
%
% Outputs:
%
%   fid            - raw column of coil'*rho values, starting just after
%                    pulse 1, kick 1, pulse 2, and kick 2 (no delays)
%
% Requires numeric imaging generators, a Liouville-space basis, and N>1
% equal cells. Kicks use the bare signed isotope gyromagnetic ratio at
% cell centres, not the site-specific shielding of a finite gradient
% lobe. RF uses exp(-1i*angle*Lx); no receiver rephasing, volume
% normalisation, or initial-point rescaling is applied. Include volume
% weights in the receive phantom. Check spatial and time-step convergence
% separately: nonlinear mixing can generate unresolved spatial harmonics.

function fid=crazed_ddf(spin_system,parameters,H,R,K,~,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F);

% Cache the spin commutators and the corresponding unweighted probes
isotope=parameters.spins{1};
labels={'Lx','Ly','Lz'};
ops=cell(1,3);
probe_cols=cell(1,3);
for n=1:3
    ops{n}=operator(spin_system,labels{n},isotope);
    probe_cols{n}=coil_state(spin_system,labels{n},isotope);
end

% Calibrate the longitudinal projection once against the declared reference
probes=[probe_cols{:}];
ref=probes(:,3)'*parameters.rho_ref;
if (~isfinite(ref))||(real(ref)==0)||...
   (abs(imag(ref))>spin_system.tols.liouv_zero*abs(real(ref)))
    error('rho_ref must have a real nonzero longitudinal projection.');
end
scale=parameters.magnetisation/real(ref);

% Check periodicity without rounding the imposed area
gamma=spin(isotope);
turns=gamma*parameters.gradient_areas*parameters.dims/(2*pi);
if any(abs(turns-round(turns))>8*eps(max(1,abs(turns))))
    error('each gradient area must encode whole pitches over the sample.');
end

% Construct the bare Zeeman ramp on the periodic cell centres
npts=parameters.npts;
coord=((0:npts-1)'-(npts-1)/2)*(parameters.dims/npts);
gradient=kron(spdiags(-gamma*coord,0,npts,npts),ops{3});
pulse=kron(speye(npts),ops{1});

% Apply the instantaneous X pulses and their signed gradient kicks
rho=parameters.rho0;
for n=1:2
    rho=step(spin_system,pulse,rho,parameters.flip_angles(n));
    rho=step(spin_system,gradient,rho,parameters.gradient_areas(n));
end

% Acquire the fixed receiver at the end of preparation
fixed=sparse(H+F+1i*R+1i*K);
fid=zeros(parameters.npoints,1);
fid(1)=parameters.coil'*rho;
dt=1/(parameters.sweep*parameters.ddf_steps);

% Propagate between samples with the field evaluated at every LG4 stage
for n=2:parameters.npoints

    % Evolve either the linear control or the self-consistent state
    if parameters.magnetisation==0
        rho=step(spin_system,fixed,rho,1/parameters.sweep);
    else
        for k=1:parameters.ddf_steps
            time=((n-2)*parameters.ddf_steps+k-1)*dt;
            rho=step(spin_system,{@generator,time,'LG4'},rho,dt);
        end
    end

    % Detect without a state-dependent receiver phase
    fid(n)=parameters.coil'*rho;

end

    % Assemble the local Zeeman feedback from the current stage state
    function L=generator(~,stage)

        % A Hermitian state has real Cartesian moments in every voxel
        moments=transpose(probes'*reshape(stage,parameters.spn_dim,npts));
        if any(max(abs(imag(moments)),[],2)>spin_system.tols.liouv_zero*...
               max(abs(real(ref)),max(abs(real(moments)),[],2)))
            error('non-real Cartesian spin projections in the DDF source.');
        end

        % Convert the Cartesian magnetisation to tesla and then rad/s
        field=ddf_field(scale*full(real(moments)),parameters.direction,parameters.demag_factor);
        L=fixed;
        for a=1:3
            L=L+kron(spdiags(-gamma*field(:,a),0,npts,npts),ops{a});
        end

    end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('a Liouville-space basis is required.');
end
if ~all(cellfun(@isnumeric,{H,R,K,F}))
    error('numeric imaging generators are required; disable polyadic storage.');
end
if (~iscell(parameters.spins))||(~isscalar(parameters.spins))||...
   (~ischar(parameters.spins{1}))||(~all(strcmp(parameters.spins{1},spin_system.comp.isotopes)))
    error('one isotope common to all spins must be selected.');
end
if (~isnumeric(parameters.npts))||(~isscalar(parameters.npts))||(~isreal(parameters.npts))||...
   (~isfinite(parameters.npts))||(parameters.npts~=fix(parameters.npts))||(parameters.npts<=1)
    error('parameters.npts must be an integer greater than one.');
end
if (~isfloat(parameters.dims))||(~isscalar(parameters.dims))||(~isreal(parameters.dims))||...
   (~isfinite(parameters.dims))||(parameters.dims<=0)
    error('parameters.dims must be a real finite positive floating-point scalar.');
end
if ~ismember(parameters.deriv{1},{'period','fourier'})
    error('the spatial derivative must be periodic.');
end
if (~isfloat(parameters.rho_ref))||(any(~isfinite(parameters.rho_ref(:))))||...
   (~isequal(size(parameters.rho_ref),[parameters.spn_dim 1]))
    error('parameters.rho_ref must be a finite local spin-state column.');
end
if (~isfloat(parameters.rho0))||(any(~isfinite(parameters.rho0(:))))||...
   (~isequal(size(parameters.rho0),[parameters.spn_dim*parameters.npts 1]))
    error('parameters.rho0 must be a finite column matching the spin-space dimension.');
end
if (~isfloat(parameters.coil))||(any(~isfinite(parameters.coil(:))))||...
   (~isequal(size(parameters.coil),[parameters.spn_dim*parameters.npts 1]))
    error('parameters.coil must be a finite column matching the spin-space dimension.');
end
if (~isfloat(parameters.direction))||(~isvector(parameters.direction))||...
   (numel(parameters.direction)~=3)||(~isreal(parameters.direction))||...
   (any(~isfinite(parameters.direction(:))))
    error('parameters.direction must be a real finite floating-point 3-element vector.');
end
if abs(norm(parameters.direction)-1)>3*eps(class(parameters.direction))
    error('direction must be a unit vector.');
end
if (~isfloat(parameters.demag_factor))||(~isscalar(parameters.demag_factor))||...
   (~isreal(parameters.demag_factor))||(parameters.demag_factor<0)||(parameters.demag_factor>1)||...
   (~isfinite(parameters.demag_factor))
    error('parameters.demag_factor must be a real scalar in the interval [0,1].');
end
if (~isfloat(parameters.magnetisation))||(~isscalar(parameters.magnetisation))||...
   (~isreal(parameters.magnetisation))||(~isfinite(parameters.magnetisation))||...
   (parameters.magnetisation<0)
    error('parameters.magnetisation must be a real finite nonnegative floating-point scalar.');
end
if (~isfloat(parameters.flip_angles))||(~isequal(size(parameters.flip_angles),[1 2]))||...
   (~isreal(parameters.flip_angles))||(any(~isfinite(parameters.flip_angles(:))))
    error('parameters.flip_angles must be a real finite floating-point 1-by-2 array.');
end
if (~isfloat(parameters.gradient_areas))||(~isequal(size(parameters.gradient_areas),[1 2]))||...
   (~isreal(parameters.gradient_areas))||(any(~isfinite(parameters.gradient_areas(:))))
    error('parameters.gradient_areas must be a real finite floating-point 1-by-2 array.');
end
if (~isfloat(parameters.sweep))||(~isscalar(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isfinite(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a real finite positive floating-point scalar.');
end
if (~isnumeric(parameters.npoints))||(~isscalar(parameters.npoints))||...
   (~isreal(parameters.npoints))||(~isfinite(parameters.npoints))||...
   (parameters.npoints~=fix(parameters.npoints))||(parameters.npoints<=0)
    error('parameters.npoints must be a positive integer.');
end
if (~isnumeric(parameters.ddf_steps))||(~isscalar(parameters.ddf_steps))||...
   (~isreal(parameters.ddf_steps))||(~isfinite(parameters.ddf_steps))||...
   (parameters.ddf_steps~=fix(parameters.ddf_steps))||(parameters.ddf_steps<=0)
    error('parameters.ddf_steps must be a positive integer.');
end
end


