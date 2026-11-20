% Piecewise-constant imaging sequence with a self-consistent secular
% distant dipolar field in separated rotating frames. Syntax:
%
%                 fid=imaging(spin_system,@ddf_sequence,p)
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

function fid=ddf_sequence(spin_system,p,H,R,K,G,F)

% Check consistency
grumble(spin_system,p,H,R,K,G,F);

% Prepare the spin-space rotations and local unweighted probes
labels={'Lx','Ly','Lz'};
ngroups=numel(p.groups);
ops=cell(ngroups,3);
probes=cell(ngroups,1);
gamma=zeros(ngroups,1);
for q=1:ngroups
    gamma(q)=spin(spin_system.comp.isotopes{p.groups{q}(1)});
    cols=cell(1,3);
    for a=1:3
        ops{q,a}=operator(spin_system,labels{a},p.groups{q});
        cols{a}=coil_state(spin_system,labels{a},p.groups{q});
    end
    probes{q}=[cols{:}];
end

% Convert endpoint-to-endpoint context ramps into cell-centred ramps
for a=1:3
    G{a}=G{a}*((p.npts(a)-1)/p.npts(a));
end

% Assemble the common drift and initialise the receiver
fixed=sparse(H+F+p.transport+1i*R+1i*K);
rho=p.rho0;
fid=zeros(nnz(p.acquire),1);
sample=0;
time=0;

% Hold the external controls constant between waveform edges
for n=1:numel(p.duration)

    % Combine group-selective RF and three simultaneous gradient axes
    segment=fixed;
    for a=1:3
        segment=segment+p.gradient(n,a)*G{a};
    end
    for q=1:ngroups
        pulse=real(p.rf(n,q))*ops{q,1}+imag(p.rf(n,q))*ops{q,2};
        segment=segment+kron(speye(p.spc_dim),pulse);
    end

    % Evaluate the field at the internal stages of each nonlinear step
    if p.duration(n)>0
        if any(p.scale)
            dt=p.duration(n)/p.steps(n);
            for k=1:p.steps(n)
                rho=step(spin_system,{@generator,time+(k-1)*dt,'LG4'},rho,dt);
            end
        else
            rho=step(spin_system,segment,rho,p.duration(n));
        end
    end
    time=time+p.duration(n);

    % Sample only the declared receiver without changing its phase
    if p.acquire(n)
        sample=sample+1;
        fid(sample)=p.coil'*rho;
    end

end

    % Assemble the stage-dependent Zeeman term in each target frame
    function L=generator(~,stage)

        % Sum the magnetisation in SI units separately for each frame
        states=reshape(stage,p.spn_dim,p.spc_dim);
        L=segment;
        for frame=unique(p.frames)

            % Project local moments without weighting by the receive coil
            magn=zeros(p.spc_dim,3);
            for b=find((p.frames==frame)&(p.scale~=0))
                moments=transpose(probes{b}'*states);
                bound=norm(probes{b},'fro')*vecnorm(states).';
                if any(abs(imag(moments))>spin_system.tols.liouv_zero*bound,'all')
                    error('non-real Cartesian magnetisation in a source group.');
                end
                magn=magn+p.scale(b)*real(moments);
            end
            if ~any(magn,'all'), continue; end

            % Couple the field into each target with its own Zeeman factor
            field=p.field(magn);
            validateattributes(field,{'double','single'},{'real','finite','size',[p.spc_dim 3]});
            for t=1:ngroups
                axes=3;
                if p.frames(t)==frame, axes=1:3; end
                for c=axes
                    L=L+kron(spdiags(-gamma(t)*field(:,c),0,p.spc_dim,p.spc_dim),ops{t,c});
                end
            end

        end

    end

end

% Consistency enforcement
function grumble(spin_system,p,H,R,K,G,F)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('a Liouville-space basis is required.');
end
if ~all(cellfun(@isnumeric,[{H,R,K,F} G]))
    error('numeric imaging generators are required.');
end
validateattributes(p.npts,{'numeric'},{'size',[1 3],'integer','positive'});
if isfield(p,'grad_angles')
    error('provide the gradient components in unrotated laboratory axes.');
end
validateattributes(p.groups,{'cell'},{'row','nonempty'});
ngroups=numel(p.groups);
validateattributes(p.frames,{'numeric'},{'size',[1 ngroups],'real','integer','finite'});
validateattributes(p.scale,{'double','single'},{'size',[1 ngroups],'real','finite'});
for q=1:ngroups
    validateattributes(p.groups{q},{'numeric'},{'row','nonempty','integer','positive','<=',spin_system.comp.nspins});
end
if ~isequal(sort([p.groups{:}]),1:spin_system.comp.nspins)
    error('groups must partition the spin sites.');
end
for frame=unique(p.frames)
    sites=[p.groups{p.frames==frame}];
    if numel(unique(spin_system.comp.isotopes(sites)))~=1
        error('a common frame must contain a single isotope.');
    end
end
validateattributes(p.field,{'function_handle'},{'scalar'});
dim=p.spn_dim*p.spc_dim;
validateattributes(p.transport,{'double'},{'size',[dim dim],'finite'});
if ~issparse(p.transport), error('transport must be sparse.'); end
validateattributes(p.rho0,{'double','single'},{'size',[dim 1],'finite'});
validateattributes(p.coil,{'double','single'},{'size',[dim 1],'finite'});
validateattributes(p.duration,{'double','single'},{'column','nonempty','real','finite','nonnegative'});
count=numel(p.duration);
validateattributes(p.rf,{'double','single'},{'size',[count ngroups],'finite'});
validateattributes(p.gradient,{'double','single'},{'size',[count 3],'real','finite'});
validateattributes(p.steps,{'numeric'},{'size',[count 1],'integer','positive','finite'});
validateattributes(p.acquire,{'logical'},{'size',[count 1]});
end


