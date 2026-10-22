% Field swept powder EPR spectra. A rough implementation with ex-
% pensive eigenfields algorithm, an explicit spherical grid, and
% a hard-coded Lorentzian line shape. Syntax:
%
%        [spec,parameters]=fieldsweep(spin_system,parameters)
%
% Parameters:
%
%     parameters.grid      -  initial spherical grid, ideally
%                             a non-symmetric one to avoid
%                             transition degeneracies, a good
%                             start is
%                             'rep_2ang_100pts_sph'
%
%     parameters.spins     -  a one-element cell array speci-
%                             fying the spin that is coupled
%                             couple to the microwave field,
%                             e.g. {'E'}
%
%     parameters.mw_freq   -  microwave frequency, Hz
%
%     parameters.fwhm      -  Lorentzian line FWHM, Tesla
%
%     parameters.window    -  field sweep window in Tesla,
%                             as a vector [Bmin Bmax]
%
%     parameters.npoints   -  number of points in the sweep
%
%     parameters.tm_tol    -  relative transition moment to-
%                             lerance, 0.01 is a good start
%
%     parameters.rspt_order - perturbation theory order for
%                             eigenfields calculation, 2 is
%                             a good start; specify Inf for
%                             exact diagonalisation
%
%     parameters.int_tol   -  powder integration tolerance,
%                             a balance between speed and
%                             integration accuracy
%
% Outputs:
%
%     spec       - field-swept EPR spectrum
%
%     parameters - updated parameters structure with the mag-
%                  netic field axis for plotting returned in
%                  parameters.b_axis
%
% Note: irrespective of the actual sweep extents, the magnetic field
%       in sys.magnet should be set to 1 Tesla.
%
% Note: this experiment should be called directly without a context.
%       Adaptive refinement is scheduled on the client. Shared edges
%       are reserved while in flight, then cached for neighbouring
%       triangles; each midpoint eigenset is computed only once.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fieldsweep.m>

function [spec,parameters]=fieldsweep(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Peak position tolerance is quarter a grid pixel
window_size=max(parameters.window)-min(parameters.window);
parameters.pp_tol=0.25*window_size/(parameters.npoints-1);

% Get the Hamiltonians and tidy up their isotropic parts
[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));
[Iz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));
Ic=(Ic+Ic')/2; Iz=(Iz+Iz')/2; 

% Get the microwave operator (state in Liouville space)
Hmw=coil_state(spin_system,'Lx',parameters.spins{1});

% Get the initial grid and compute its convex hull
load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' ... 
      filesep parameters.grid],'betas','gammas');
hull=get_hull(betas,gammas); grid_size=numel(betas);

% Make the magnetic field axis
parameters.b_axis=linspace(parameters.window(1),...
                           parameters.window(2),...
                           parameters.npoints);

% Transfer the heavy operators once and reuse them in every task
operators=parallel.pool.Constant({Ic,Iz,Qc,Qz,Hmw});

% Schedule the initial vertex eigensets with worker-local operators
futures(1:grid_size)=parallel.Future;
for n=1:grid_size
    localpar=parameters; localpar.orientation=[0 betas(n) gammas(n)];
    futures(n)=parfeval(@vertex_task,1,spin_system,localpar,operators);
end

% Retrieve vertices and attach their spherical coordinates
eigensets=fetchOutputs(futures);
for n=1:grid_size
    eigensets(n).xyz=[sin(betas(n))*cos(gammas(n));...
                      sin(betas(n))*sin(gammas(n));...
                      cos(betas(n))];
end

% Index cached nodes by sorted edge endpoints; empty nodes are in flight
nodes=num2cell(eigensets); triangles=hull;
edges=containers.Map('KeyType','char','ValueType','double');
midpoint_ids=zeros(size(triangles));
futures=parallel.FevalFuture.empty(1,0); active=zeros(1,0);
pending=(1:size(triangles,1))';
spec=zeros(size(parameters.b_axis),'like',1i);

% Refine on the client, allowing branches to use the whole pool
while (~isempty(pending))||(~isempty(futures))

    % Submit triangles whose shared edges are not already in flight
    waiting=false(size(pending));
    for p=1:numel(pending)
        n=pending(p);

        % Look up the three edges in the midpoint cache
        vertices=triangles(n,:); pairs=vertices([1 2; 2 3; 3 1]);
        keys=cell(1,3); ids=zeros(1,3); blocked=false;
        for k=1:3
            keys{k}=sprintf('%d,%d',sort(pairs(k,:)));
            if isKey(edges,keys{k})
                ids(k)=edges(keys{k});
                blocked=blocked||isempty(nodes{ids(k)});
            end
        end

        % A reserved edge must finish before another task can use it
        if blocked
            waiting(p)=true; continue
        end

        % Reserve missing midpoints before another triangle is submitted
        midpoints=cell(1,3);
        for k=1:3
            if ids(k)==0
                nodes{end+1}=[]; ids(k)=numel(nodes);
                edges(keys{k})=ids(k);
            else
                midpoints{k}=nodes{ids(k)};
            end
        end
        midpoint_ids(n,:)=ids;

        % Evaluate one refinement level using cached midpoint eigensets
        triangle=[nodes{vertices}];
        futures(end+1)=parfeval(@triangle_task,3,spin_system,...
                               parameters,triangle,operators,midpoints);
        active(end+1)=n;

    end
    pending=pending(waiting);

    % Drain the next completed triangle and release its reserved edges
    [idx,local_spec,err,midpoints]=fetchNext(futures);
    n=active(idx); ids=midpoint_ids(n,:);
    nodes(ids)=midpoints;
    futures(idx)=[]; active(idx)=[];

    % Accumulate converged triangles or expose four independent children
    if err>parameters.int_tol
        vertices=triangles(n,:); first=size(triangles,1)+1;
        children=[vertices(1) ids(1) ids(3);...
                  ids(1) vertices(2) ids(2);...
                  ids(3) ids(2) vertices(3);...
                  ids(1) ids(2) ids(3)];
        triangles=[triangles; children];
        midpoint_ids=[midpoint_ids; zeros(4,3)];
        pending=[(first:first+3)'; pending];
    else
        spec=spec+local_spec;
    end

end

end

% Assemble an initial vertex eigenset from worker-local operators
function eigenset=vertex_task(spin_system,parameters,operators)

% Orient the operators at the requested vertex
operators=operators.Value;
Hz=operators{2}+orientation(operators{4},parameters.orientation);
Hc=operators{1}+orientation(operators{3},parameters.orientation);

% Solve the field-swept eigenproblem
eigenset=eigenfields(spin_system,parameters,Hz,Hc,operators{5});

end

% Evaluate a triangle using the operators resident on this worker
function [spec,err,midpoints]=triangle_task(spin_system,parameters,...
                                           triangle,operators,midpoints)

% Unwrap the constant without capturing the operators in a closure
operators=operators.Value;
[spec,err,midpoints]=voitlander(spin_system,parameters,triangle,...
                              operators{:},midpoints);

end

% Consistency enforcement
function grumble(spin_system,parameters)
if spin_system.inter.magnet~=1
    error('field swept experiment, set sys.magnet=1');
end
if ~isfield(parameters,'grid')
    error('grid must be specified in parameters.grid variable.');
end
if ~ischar(parameters.grid)
    error('parameters.grid must be character string.');
end
if ~isfield(parameters,'spins')
    error('spin to couple to the microwave field must be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
   (~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
end
if ~isfield(parameters,'mw_freq')
    error('MW frequency should be specified in parameters.mw_freq variable.');
end
if (~isnumeric(parameters.mw_freq))||(~isscalar(parameters.mw_freq))||...
   (~isreal(parameters.mw_freq))
    error('parameters.mw_freq must be a real scalar.');
end
if ~isfield(parameters,'fwhm')
    error('line width should be specified in parameters.fwhm variable.');
end
if (~isnumeric(parameters.fwhm))||(~isscalar(parameters.fwhm))||...
   (~isreal(parameters.fwhm))||(parameters.fwhm<=0)
    error('parameters.fwhm must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isscalar(parameters.npoints))||...
   (~isreal(parameters.npoints))||(parameters.npoints<=1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a real integer greater than one.');
end
if ~isfield(parameters,'window')
    error('sweep window should be specified in parameters.window variable.');
end
if (~isnumeric(parameters.window))||(numel(parameters.window)~=2)||...
   (~isreal(parameters.window))||(parameters.window(2)<=parameters.window(1))
    error('parameters.window must contain two ascending numbers.');
end
if ~isfield(parameters,'tm_tol')
    error('transition moment tolerance must be specified in parameters.tm_tol variable.');
end
if (~isnumeric(parameters.tm_tol))||(~isscalar(parameters.tm_tol))||...
   (~isreal(parameters.tm_tol))||(parameters.tm_tol<0)
    error('parameters.tm_tol must be a non-negative real scalar.');
end
if ~isfield(parameters,'rspt_order')
    error('perturbation theory order must be specified in parameters.rspt_order variable.');
end
if (~isnumeric(parameters.rspt_order))||(~isscalar(parameters.rspt_order))||...
   (~isreal(parameters.rspt_order))||((mod(parameters.rspt_order,1)~=0)&&...
   (~isinf(parameters.rspt_order)))||(parameters.rspt_order<0)
    error('parameters.rspt_order must be a non-negative integer or Inf.');
end
if ~isfield(parameters,'int_tol')
    error('integration tolerance must be specified in parameters.int_tol variable.');
end
if (~isnumeric(parameters.int_tol))||(~isscalar(parameters.int_tol))||...
   (~isreal(parameters.int_tol))||(parameters.int_tol<=0)
    error('parameters.int_tol must be a positive real scalar.');
end
end

% He leads the existence of a real bohemian intellectual. Washing, 
% grooming, and changing his linen are things he does rarely, and
% he likes to get drunk. Though he is often idle for days on end,
% he will work day and night with tireless endurance when he has 
% a great deal of work to do. He has no fixed times for going to 
% sleep and waking up. He often stays up all night, and then lies
% down fully clothed on the sofa at midday and sleeps till evening,
% untroubled by the comings and goings of the whole world.
%
% From the Prussian secret police
% report on Karl Marx


