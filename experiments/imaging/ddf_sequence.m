% Piecewise-constant imaging sequence with a self-consistent secular
% distant dipolar field in separated rotating frames. Syntax:
%
%                 fid=imaging(spin_system,@ddf_sequence,parameters)
%
% Parameters (in addition to ordinary numeric 3D imaging inputs):
%
%   groups     - 1-by-Q cell of site-index rows partitioning the spins;
%                each group must contain only one isotope
%   frames     - 1-by-Q integer labels; groups in a common rotating
%                frame share transverse fields and must be one isotope
%   scale      - 1-by-Q signed A/m per unit Cartesian moment projection;
%                include population only if not already in the state;
%                zero makes a group a passive target
%   field      - linear field handle mapping N-by-3 real magnetisation
%                densities (A/m) to N-by-3 real secular fields (T),
%                with voxels in MATLAB [X Y Z] linear order
%   transport  - additional sparse spin-space generator (rad/s); for
%                independent molecules use 1i times their spatial rate
%                matrices lifted onto the respective substance blocks;
%                set context diffusion/flow to zero when replacing F
%   duration   - S-by-1 segment durations (s); zero allows initial read
%   gradient   - S-by-3 Cartesian gradients (T/m) on voxel centres
%   rf         - S-by-Q nutation frequencies (rad/s), real X, imag Y
%   steps      - S-by-1 positive integer substep counts
%   acquire    - S-by-1 logical enabling detection at segment ends
%
% Outputs:
%
%   fid        - complex column of raw coil'*rho samples in time order
%
% Longitudinal fields act on all targets with their signed gyromagnetic
% ratios. Transverse fields act only within the source rotating frame.
% Controls, transport, and feedback act simultaneously during each
% finite segment; split the waveform wherever any control changes.

function fid=ddf_sequence(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Prepare the spin-space rotations and local unweighted probes
labels={'Lx','Ly','Lz'};
ngroups=numel(parameters.groups);
ops=cell(ngroups,3);
probes=cell(ngroups,1);
gamma=zeros(ngroups,1);
for q=1:ngroups
    gamma(q)=spin(spin_system.comp.isotopes{parameters.groups{q}(1)});
    cols=cell(1,3);
    for a=1:3
        ops{q,a}=operator(spin_system,labels{a},parameters.groups{q});
        cols{a}=coil_state(spin_system,labels{a},parameters.groups{q});
    end
    probes{q}=[cols{:}];
end

% Convert endpoint-to-endpoint context ramps into cell-centred ramps
for a=1:3
    G{a}=G{a}*((parameters.npts(a)-1)/parameters.npts(a));
end

% Assemble the common drift and initialise the receiver
fixed=sparse(H+F+parameters.transport+1i*R+1i*K);
rho=parameters.rho0;
fid=zeros(nnz(parameters.acquire),1);
sample=0;
time=0;

% Hold the external controls constant between waveform edges
for n=1:numel(parameters.duration)

    % Combine group-selective RF and three simultaneous gradient axes
    segment=fixed;
    for a=1:3
        segment=segment+parameters.gradient(n,a)*G{a};
    end
    for q=1:ngroups
        pulse=real(parameters.rf(n,q))*ops{q,1}+imag(parameters.rf(n,q))*ops{q,2};
        segment=segment+kron(speye(parameters.spc_dim),pulse);
    end

    % Evaluate the field at the internal stages of each nonlinear step
    if parameters.duration(n)>0
        if any(parameters.scale)
            dt=parameters.duration(n)/parameters.steps(n);
            for k=1:parameters.steps(n)
                rho=step(spin_system,{@generator,time+(k-1)*dt,'LG4'},rho,dt);
            end
        else
            rho=step(spin_system,segment,rho,parameters.duration(n));
        end
    end
    time=time+parameters.duration(n);

    % Sample only the declared receiver without changing its phase
    if parameters.acquire(n)
        sample=sample+1;
        fid(sample)=parameters.coil'*rho;
    end

end

    % Assemble the stage-dependent Zeeman term in each target frame
    function L=generator(~,stage)

        % Sum the magnetisation in SI units separately for each frame
        states=reshape(stage,parameters.spn_dim,parameters.spc_dim);
        L=segment;
        for frame=unique(parameters.frames)

            % Project local moments without weighting by the receive coil
            magn=zeros(parameters.spc_dim,3);
            for b=find((parameters.frames==frame)&(parameters.scale~=0))
                moments=transpose(probes{b}'*states);
                bound=norm(probes{b},'fro')*vecnorm(states).';
                if any(abs(imag(moments))>spin_system.tols.liouv_zero*bound,'all')
                    error('non-real Cartesian magnetisation in a source group.');
                end
                magn=magn+parameters.scale(b)*real(moments);
            end
            if ~any(magn,'all')
                continue
            end

            % Couple the field into each target with its own Zeeman factor
            field=parameters.field(magn);
            if (~isfloat(field))||(~isreal(field))||(any(~isfinite(field(:))))||...
               (~isequal(size(field),[parameters.spc_dim 3]))
                error('field must be a finite real array of three Cartesian components per voxel.');
            end
            for t=1:ngroups
                axes=3;
                if parameters.frames(t)==frame
                    axes=1:3;
                end
                for c=axes
                    L=L+kron(spdiags(-gamma(t)*field(:,c),0,parameters.spc_dim,parameters.spc_dim),ops{t,c});
                end
            end

        end

    end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('a Liouville-space basis is required.');
end
if ~all(cellfun(@isnumeric,[{H,R,K,F} G]))
    error('numeric imaging generators are required.');
end
if (~isnumeric(parameters.npts))||(~isequal(size(parameters.npts),[1 3]))||...
   (~isreal(parameters.npts))||(any(~isfinite(parameters.npts(:))))||...
   (any(parameters.npts(:)~=fix(parameters.npts(:))))||(any(parameters.npts(:)<=0))
    error('parameters.npts must be a 1-by-3 row of positive integers.');
end
if isfield(parameters,'grad_angles')
    error('provide the gradient components in unrotated laboratory axes.');
end
if (~iscell(parameters.groups))||(~isrow(parameters.groups))||(isempty(parameters.groups))
    error('parameters.groups must be a non-empty cell row.');
end
ngroups=numel(parameters.groups);
if (~isnumeric(parameters.frames))||(~isequal(size(parameters.frames),[1 ngroups]))||...
   (~isreal(parameters.frames))||(any(~isfinite(parameters.frames(:))))||...
   (any(parameters.frames(:)~=fix(parameters.frames(:))))
    error('parameters.frames must be a row containing one finite integer frame label per group.');
end
if (~isfloat(parameters.scale))||(~isequal(size(parameters.scale),[1 ngroups]))||...
   (~isreal(parameters.scale))||(any(~isfinite(parameters.scale(:))))
    error('parameters.scale must be a row containing one finite real magnetisation scale per group.');
end
for q=1:ngroups
    if (~isnumeric(parameters.groups{q}))||(~isrow(parameters.groups{q}))||...
       (isempty(parameters.groups{q}))||(~isreal(parameters.groups{q}))||...
       (any(~isfinite(parameters.groups{q}(:))))||...
       (any(parameters.groups{q}(:)~=fix(parameters.groups{q}(:))))||...
       (any(parameters.groups{q}(:)<=0))||(any(parameters.groups{q}(:)>spin_system.comp.nspins))
        error('parameters.groups{q} must be a non-empty row of spin indices between 1 and the number of spins.');
    end
end
if ~isequal(sort([parameters.groups{:}]),1:spin_system.comp.nspins)
    error('groups must partition the spin sites.');
end
for frame=unique(parameters.frames)
    sites=[parameters.groups{parameters.frames==frame}];
    if numel(unique(spin_system.comp.isotopes(sites)))~=1
        error('a common frame must contain a single isotope.');
    end
end
if (~isa(parameters.field,'function_handle'))||(~isscalar(parameters.field))
    error('parameters.field must be a scalar function handle.');
end
dim=parameters.spn_dim*parameters.spc_dim;
if (~isa(parameters.transport,'double'))||(~isequal(size(parameters.transport),[dim dim]))||...
   (any(~isfinite(parameters.transport(:))))
    error('parameters.transport must be a finite double matrix on the full spin-space basis.');
end
if ~issparse(parameters.transport)
    error('parameters.transport must be sparse.');
end
if (~isfloat(parameters.rho0))||(~isequal(size(parameters.rho0),[dim 1]))||...
   (any(~isfinite(parameters.rho0(:))))
    error('parameters.rho0 must be a finite column matching the spin-space dimension.');
end
if (~isfloat(parameters.coil))||(~isequal(size(parameters.coil),[dim 1]))||...
   (any(~isfinite(parameters.coil(:))))
    error('parameters.coil must be a finite column matching the spin-space dimension.');
end
if (~isfloat(parameters.duration))||(~iscolumn(parameters.duration))||...
   (isempty(parameters.duration))||(~isreal(parameters.duration))||...
   (any(~isfinite(parameters.duration(:))))||(any(parameters.duration(:)<0))
    error('parameters.duration must be a non-empty real finite nonnegative floating-point column.');
end
count=numel(parameters.duration);
if (~isfloat(parameters.rf))||(~isequal(size(parameters.rf),[count ngroups]))||...
   (any(~isfinite(parameters.rf(:))))
    error('parameters.rf must be a finite RF array with one row per segment and one column per group.');
end
if (~isfloat(parameters.gradient))||(~isequal(size(parameters.gradient),[count 3]))||...
   (~isreal(parameters.gradient))||(any(~isfinite(parameters.gradient(:))))
    error('parameters.gradient must be a finite real array of three gradient components per segment.');
end
if (~isnumeric(parameters.steps))||(~isequal(size(parameters.steps),[count 1]))||...
   (~isreal(parameters.steps))||(any(~isfinite(parameters.steps(:))))||...
   (any(parameters.steps(:)~=fix(parameters.steps(:))))||(any(parameters.steps(:)<=0))
    error('parameters.steps must be a column of positive substep counts, one per waveform segment.');
end
if (~islogical(parameters.acquire))||(~isequal(size(parameters.acquire),[count 1]))
    error('parameters.acquire must be a logical column with one entry per segment.');
end
end


