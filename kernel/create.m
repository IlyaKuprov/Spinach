% The entry function of the Spinach kernel that creates the spin 
% system object that the rest of the library requires to run. It
% checks and absorbs interaction specifications, and writes some
% useful diagnostics to the console. Syntax:
% 
%                  spin_system=create(sys,inter)
%
% Parameters:
%
%    sys   - spin system and instrument specification
%            structure, see the spin system specifica-
%            tion section of the online manual
%
%  inter   - interaction specification structure, see
%            see the spin system specification section
%            of the online manual
%
% Outputs:
%
%  spin_system   - the primary object used by Spinach
%                  to store simulation information
%
% i.kuprov@soton.ac.uk
% hannah.hogben@chem.ox.ac.uk
% kpervushin@ntu.edu.sg
% luke.edwards@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=create.m>

function spin_system=create(sys,inter)

% Refuse to run in script mode
call_stack=dbstack;
if strcmp(call_stack(end).name,'create')&&(~isworkernode)
    error('Never run Matlab in script mode. All your variables persist, and\n%s',...
          'there will be no end to confusion. Add a [function ... end] block.');
end

% Close all open files
fclose('all');

% Locate the root and run sanity checks
if isempty(which('existentials'))
    
    % Tell the user to RTFM
    error('paths have not been correctly set - please follow installation\n%s',...
          'instructions (make sure you had included the subdirectories).');
       
else
    
    % Set the root directory
    root_dir=which('existentials');
    spin_system.sys.root_dir=root_dir(1:(end-32));
    
end

% Run overrides
autoexec;

% Validate input
grumble(sys,inter);

% Run existential checks
if (~isfield(sys,'disable'))||...
   (~ismember('hygiene',sys.disable))
    existentials();
end

% Decide output destination
if isfield(sys,'output')
    
    % String specifications
    if strcmp(sys.output,'hush')
    
        % Hush the output
        spin_system.sys.output='hush';
        
    elseif strcmp(sys.output,'console')
    
        % Print to the console
        spin_system.sys.output=1;
        
    else
    
        % Print to a user-specified file
        spin_system.sys.output=fopen(sys.output,'a');
        
    end
    
    % Parse out
    sys=rmfield(sys,'output');
    
else
    
    % Print to the console
    spin_system.sys.output=1;
        
end

% Decide scratch destination
if isfield(sys,'scratch')
    
    % Scratch to a user-specified directory
    spin_system.sys.scratch=sys.scratch;
    
    % Parse out
    sys=rmfield(sys,'scratch');
    
else
    
    % Scratch to the default directory
    spin_system.sys.scratch=[spin_system.sys.root_dir filesep 'scratch'];
    
end

% If scratch directory does not exist, create it
if ~exist(spin_system.sys.scratch,'dir')
    mkdir(spin_system.sys.scratch);
end

% Report scratch location
report(spin_system,['scratch location: ' spin_system.sys.scratch]);

% Show version banner
banner(spin_system,'version_banner');

% Internal algorithms to disable
if isfield(sys,'disable')
    
    % Absorb user-specified values
    spin_system.sys.disable=sys.disable;
    
    % Parse out
    sys=rmfield(sys,'disable');
    
else
    
    % Disable nothing by default
    spin_system.sys.disable={};
    
end

% Internal algorithms to enable
if isfield(sys,'enable')
    
    % Absorb user-specified values
    spin_system.sys.enable=sys.enable;
    
    % Parse out
    sys=rmfield(sys,'enable');
    
else
    
    % Enable nothing by default
    spin_system.sys.enable={};
    
end

% Process internal tolerances
[spin_system,sys]=tolerances(spin_system,sys);

% Tidy up the cache
if ~isworkernode, cacheman(spin_system); end

% Disabled features report
if ~isempty(spin_system.sys.disable)
    report(spin_system,'WARNING: the following functionality is disabled by the user');
    if ismember('hygiene',spin_system.sys.disable),   report(spin_system,'         > health checks at start-up'); end
    if ismember('pt',spin_system.sys.disable),        report(spin_system,'         > automatic detection of non-interacting subspaces'); end
    if ismember('zte',spin_system.sys.disable),       report(spin_system,'         > automatic elimination of unpopulated states'); end
    if ismember('symmetry',spin_system.sys.disable),  report(spin_system,'         > permutation symmetry factorisation'); end
    if ismember('krylov',spin_system.sys.disable),    report(spin_system,'         > Krylov propagation inside evolution() function'); end
    if ismember('clean-up',spin_system.sys.disable),  report(spin_system,'         > sparse array clean-up'); end
    if ismember('dss',spin_system.sys.disable),       report(spin_system,'         > destination state screening inside evolution() function'); end
    if ismember('expv',spin_system.sys.disable),      report(spin_system,'         > Krylov propagation inside step() function'); end
    if ismember('trajlevel',spin_system.sys.disable), report(spin_system,'         > trajectory analysis inside evolution() function'); end
    if ismember('merge',spin_system.sys.disable),     report(spin_system,'         > small subspace merging inside evolution() function'); end
    if ismember('colorbar',spin_system.sys.disable),  report(spin_system,'         > colorbar drawing by plotting utilities'); end
    if ismember('asyredf',spin_system.sys.disable),   report(spin_system,'         > asynchronous Redfied superoperator evaluation'); end
end

% Enabled features report
if ~isempty(spin_system.sys.enable)
    report(spin_system,'WARNING: the following functionality is enabled by the user');
    if ismember('gpu',spin_system.sys.enable),        report(spin_system,'         > GPU arithmetic'); end
    if ismember('op_cache',spin_system.sys.enable),   report(spin_system,'         > operator caching'); end
    if ismember('prop_cache',spin_system.sys.enable), report(spin_system,'         > propagator caching'); end
    if ismember('greedy',spin_system.sys.enable),     report(spin_system,'         > greedy parallelisation'); end
    if ismember('xmemlist',spin_system.sys.enable),   report(spin_system,'         > state-cluster cross-membership list generation'); end
    if ismember('paranoia',spin_system.sys.enable),   report(spin_system,'         > paranoid numerical accuracy settings'); end
    if ismember('cowboy',spin_system.sys.enable),     report(spin_system,'         > loose numerical accuracy settings'); end
    if ismember('polyadic',spin_system.sys.enable),   report(spin_system,'         > polyadic arithmetic with spatial degrees of freedom'); end
    if ismember('dafuq',spin_system.sys.enable),      report(spin_system,'         > detailed parallel profiling'); end
end

% Get a unique job identifier
spin_system.sys.job_id=md5_hash([clock feature('getpid')]);      %#ok<CLOCK> 
report(spin_system,['job identifier: ' spin_system.sys.job_id]);

% Create a unique scratch subdirectory for the current job
spin_system.sys.job_dir=[spin_system.sys.scratch filesep 'spinach_job_' ...
                         spin_system.sys.job_id];
if ~isfolder(spin_system.sys.job_dir), mkdir(spin_system.sys.job_dir); end

% Set up head node
if ~isworkernode
    
    % Shuffle the random number generator
    rng('shuffle'); report(spin_system,'random number generator shuffled');
    
    % Detailed parallel profiling
    if ismember('dafuq',spin_system.sys.enable)
        
        % Dump MDCS messages to console
        setSchedulerMessageHandler(@disp);
        
        % Turn on MDCS debugging
        setenv('MDCE_DEBUG','true');
        
        % Preserve job information
        pctconfig('preservejobs',true);
        
    end
    
    % Set default parallel profile
    if ~isfield(sys,'parallel')

        % Leave one core to the operating system
        sys.parallel={'local',max([1, feature('numcores')-1])};

    end
    
    % Set default parallel properties
    if ~isfield(sys,'parprops'), sys.parprops={}; end
    
    % Set up parallel pool
    current_pool=gcp('nocreate');
    if ~isempty(current_pool)
        
        % Use the existing pool
        spin_system.sys.parallel{1}=current_pool.Cluster.Profile;
        spin_system.sys.parallel{2}=current_pool.NumWorkers;
        report(spin_system,'a running parallel pool found:')
        report(spin_system,['         > parallel profile: ' spin_system.sys.parallel{1}]);
        report(spin_system,['         > worker processes: ' num2str(spin_system.sys.parallel{2})]);
        report(spin_system,'WARNING - this pool will be re-used');
        sys=rmfield(sys,'parallel'); sys=rmfield(sys,'parprops');

        % Empty ValueStore of the current pool
        store=current_pool.ValueStore; store.remove(store.keys);
        report(spin_system,'parallel pool ValueStore cleared.');

        % Disable pool timeout
        current_pool.IdleTimeout=inf;
        
    else

        % Absorb parallel profile
        spin_system.sys.parallel=sys.parallel; sys=rmfield(sys,'parallel');
        spin_system.sys.parprops=sys.parprops; sys=rmfield(sys,'parprops');
        report(spin_system,'starting a parallel pool with the following parameters:')
        report(spin_system,['         > parallel profile: ' spin_system.sys.parallel{1}]);
        report(spin_system,['         > workers to start: ' num2str(spin_system.sys.parallel{2})]);

        % Get cluster object
        c=parcluster(spin_system.sys.parallel{1});
        
        % Point the cluster into the job directory
        c.JobStorageLocation=spin_system.sys.job_dir;
        
        % Set additional properties
        for n=1:numel(spin_system.sys.parprops)
            c.AdditionalProperties.(spin_system.sys.parprops{n}{1})=...
                                    spin_system.sys.parprops{n}{2};
            if ischar(spin_system.sys.parprops{n}{2})
                report(spin_system,[' AdditionalProperties.' spin_system.sys.parprops{n}{1}...
                                    '=' spin_system.sys.parprops{n}{2}]);
            elseif isnumeric(spin_system.sys.parprops{n}{2})
                report(spin_system,[' AdditionalProperties.' spin_system.sys.parprops{n}{1}...
                                    '=' num2str(spin_system.sys.parprops{n}{2})]);
            end
        end
        
        % Start the parallel pool
        current_pool=parpool(c,spin_system.sys.parallel{2});

        % Disable pool timeout
        current_pool.IdleTimeout=inf;
    
    end
    
    % Set up parallel resource allocation
    warning('off','MATLAB:maxNumCompThreads:Deprecated');
    if ismember('greedy',spin_system.sys.enable)
        
        % Greedy
        spmd

            % Allow every worker to use every CPU core
            warning('off','MATLAB:maxNumCompThreads:Deprecated');
            maxNumCompThreads(feature('numcores'));

        end

        % Head process can use what it wants
        maxNumCompThreads(feature('numcores'));

    else
        
        % Modest
        spmd
            
            % Try to work out a reasonable allocation
            warning('off','MATLAB:maxNumCompThreads:Deprecated');
            ncores=max([1 floor(feature('numcores')/spmdSize)]);
            maxNumCompThreads(ncores);

            % First worker is special
            if spmdIndex==1
                maxNumCompThreads(feature('numcores'));
            end

        end

        % Head process can use what it wants
        maxNumCompThreads(feature('numcores'));

    end
    
end

% Make sure GPUs are present
if ismember('gpu',spin_system.sys.enable)&&(gpuDeviceCount==0)
    report(spin_system,'no CUDA capable GPUs found - running on CPU...');
    spin_system.sys.enable=setdiff(spin_system.sys.enable,{'gpu'});
end

% Binding of GPU devices to workers
if (~isworkernode)&&ismember('gpu',spin_system.sys.enable)...
                  &&(~ismember('hygiene',spin_system.sys.disable))

    % Reset all GPUs
    report(spin_system,'Clearing GPU(s)...'); 
    gpuDevice([]); pctRunOnAll('gpuDevice([]);');

    % Inform the user
    report(spin_system,'GPU configuration report:');

    % Find out how much GPU memory we have
    gpu_ram=nan(gpuDeviceCount,1);
    for n=1:gpuDeviceCount
        nth_gpu=gpuDevice(n); gpu_ram(n)=nth_gpu.TotalMemory;
        report(spin_system,['         > GPU ' num2str(n) ': ' nth_gpu.Name ...
                            ', ' num2str(gpu_ram(n)/2^30) ' GB total.']); 
    end
    
    % Head node uses the GPU with the most memory
    [~,has_most_ram]=max(gpu_ram); G=gpuDevice(has_most_ram); G.CachePolicy='minimum';
    report(spin_system,['         > head node uses GPU ' num2str(has_most_ram)]);

    % Default allocations proportional to GPU memory
    if ~isfield(sys,'gpu_bind')
        allocations=ceil(spin_system.sys.parallel{2}*gpu_ram/sum(gpu_ram));
        if sum(allocations)>spin_system.sys.parallel{2}
            allocations(has_most_ram)=allocations(has_most_ram)-1;
        end
        spin_system.sys.gpu_bind=allocations;
    else
        % Allocate manually if the user insists
        spin_system.sys.gpu_bind=sys.gpu_bind;
        sys=rmfield(sys,'gpu_bind');
    end

    % Report GPU allocations
    for n=1:gpuDeviceCount
        report(spin_system,['         > GPU ' num2str(n) ': ' ...
                            num2str(spin_system.sys.gpu_bind(n)) ' workers.']); 
    end

    % Assign workers to GPUs
    bindings=[]; report(spin_system,'Binding GPUs to workers...');
    for n=1:numel(spin_system.sys.gpu_bind)
        bindings=[bindings; n*ones(spin_system.sys.gpu_bind(n),1)]; %#ok<AGROW> 
    end
    spmd, G=gpuDevice(bindings(spmdIndex)); G.CachePolicy='minimum'; end
    
end

% Spin system banner
banner(spin_system,'spin_system_banner');

% Number and types of spins
spin_system.comp.isotopes=sys.isotopes;
spin_system.comp.nspins=numel(spin_system.comp.isotopes);
sys=rmfield(sys,'isotopes');

% Hash isotopes array for caching operations later
if ismember('op_cache',spin_system.sys.enable)
    spin_system.comp.iso_hash=md5_hash(spin_system.comp.isotopes);
end

% Text labels for spins
if isfield(sys,'labels')
    
    % Get labels from the user
    spin_system.comp.labels=sys.labels;
    sys=rmfield(sys,'labels');
    
else
    
    % Set labels to empty
    spin_system.comp.labels=cell(spin_system.comp.nspins,1);
    
end

% Multiplicities and magnetogyric ratios
spin_system.comp.mults=zeros(1,spin_system.comp.nspins);
spin_system.inter.gammas=zeros(1,spin_system.comp.nspins);
for n=1:spin_system.comp.nspins
    [spin_system.inter.gammas(n),spin_system.comp.mults(n)]=...
     spin(spin_system.comp.isotopes{n});
end
report(spin_system,['a total of ' num2str(spin_system.comp.nspins) ...
                    ' particles in the simulation, of which '...
                    num2str(nnz(spin_system.comp.mults==1)) ' have a zero spin.']);

% Primary magnet
spin_system.inter.magnet=sys.magnet; sys=rmfield(sys,'magnet');
report(spin_system,['magnetic induction of ' num2str(spin_system.inter.magnet,'%0.5g') ' Tesla ('...
                    num2str(-1e-6*spin_system.inter.magnet*spin('1H')/(2*pi),'%0.5g') ' MHz proton frequency, '...
                    num2str(-1e-9*spin_system.inter.magnet*spin('E' )/(2*pi),'%0.5g') ' GHz electron frequency).']);

% Compute carrier frequencies
spin_system.inter.basefrqs=-spin_system.inter.gammas*spin_system.inter.magnet;

% Preallocate Zeeman tensor array
spin_system.inter.zeeman.matrix=mat2cell(zeros(3*spin_system.comp.nspins,3),3*ones(spin_system.comp.nspins,1));
spin_system.inter.zeeman.ddscal=mat2cell(zeros(3*spin_system.comp.nspins,3),3*ones(spin_system.comp.nspins,1));

% Process paramagnetic shifts
if isfield(inter,'suscept')
    
    % Loop over the susceptibility centres
    for s=1:numel(inter.suscept.chi)
        
        % Get the reporting going
        report(spin_system,' ');
        report(spin_system,['nuclear paramagnetic shift summary, susceptibility centre ' num2str(s)]);
        report(spin_system,'======================================================================');
        report(spin_system,' Spin   Isotope              PMS tensor / ppm               PCS / ppm ');
        report(spin_system,'======================================================================');
    
        % Loop over the spins
        for n=1:spin_system.comp.nspins
        
            % Only process nuclei with coordinates
            if isnucleus(spin_system.comp.isotopes{n})&&(~isempty(inter.coordinates{n}))
        
                % Compute the paramagnetic shielding
                pms_tensor=xyz2pms(inter.coordinates{n},inter.suscept.xyz{s},inter.suscept.chi{s});
                
                % Report to the user
                report(spin_system,['                    ' pad(num2str(pms_tensor(1,1),'%+4.3e'),12)...
                                    pad(num2str(pms_tensor(1,2),'%+4.3e'),12) pad(num2str(pms_tensor(1,3),'%+4.3e'),12)]);
                report(spin_system,[' ' pad(num2str(n),7) pad(spin_system.comp.isotopes{n},12)...
                                    pad(num2str(pms_tensor(2,1),'%+4.3e'),12) pad(num2str(pms_tensor(2,2),'%+4.3e'),12)...
                                    pad(num2str(pms_tensor(2,3),'%+4.3e'),15) pad(num2str(trace(pms_tensor)/3,'%+4.3e'),14) ...
                                    spin_system.comp.labels{n}]);
                report(spin_system,['                    ' pad(num2str(pms_tensor(3,1),'%+4.3e'),12)...
                                    pad(num2str(pms_tensor(3,2),'%+4.3e'),12) pad(num2str(pms_tensor(3,3),'%+4.3e'),12)]);
                if n<spin_system.comp.nspins, report(spin_system,' '); end
                
                % Update the Zeeman tensor
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+pms_tensor;
                
            end
            
        end
        
        % Finish up the reporting
        report(spin_system,'======================================================================');
        
    end
    
end

% Process Zeeman interactions
if isfield(inter,'zeeman')

    % Absorb eigenvalues and Euler angles
    if isfield(inter.zeeman,'eigs')
        for n=1:spin_system.comp.nspins
            if max(abs(inter.zeeman.eigs{n}))>0
                if (~isfield(inter.zeeman,'euler'))||isempty(inter.zeeman.euler{n})
                    S=eye(3,3);
                else
                    S=euler2dcm(inter.zeeman.euler{n});
                end
                component=S*diag(inter.zeeman.eigs{n})*S';
                component=(component+component')/2;
                spin_system.inter.zeeman.matrix{n}=...
                spin_system.inter.zeeman.matrix{n}+component;
            end
        end
    end
    
    % Absorb tensors
    if isfield(inter.zeeman,'matrix')
        for n=1:spin_system.comp.nspins
            if norm(inter.zeeman.matrix{n},2)>0
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+inter.zeeman.matrix{n};
            end
        end
    end
    
    % Absorb scalars
    if isfield(inter.zeeman,'scalar')
         for n=1:spin_system.comp.nspins
            if abs(inter.zeeman.scalar{n})>0
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+eye(3)*inter.zeeman.scalar{n};
            end
        end
    end

    % Report back to the user
    summary(spin_system,'zeeman','summary of Zeeman tensors and their anisotropies (ppm for nuclei, g-tensor for electrons)');

end

% Convert Zeeman tensors to angular frequencies
if isfield(inter,'zeeman')||isfield(inter,'suscept')
    
    % Loop over the spins
    for n=1:spin_system.comp.nspins
        
        % Look at the type of the spin
        switch spin_system.comp.isotopes{n}(1)
            
            case 'E'
                
                % Store the scaling multiplier for dipolar couplings
                spin_system.inter.zeeman.ddscal{n}=spin_system.inter.zeeman.matrix{n}/spin_system.tols.freeg;
                
                % For electrons, assume that the g-tensor is given in Bohr magneton units
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}*spin_system.inter.basefrqs(n)/spin_system.tols.freeg;
                
            otherwise
                
                % Store the scaling multiplier for dipolar couplings
                spin_system.inter.zeeman.ddscal{n}=(eye(3)+(1e-6)*spin_system.inter.zeeman.matrix{n});
                
                % For nuclei, assume that the chemical shift is given in ppm
                spin_system.inter.zeeman.matrix{n}=(eye(3)+(1e-6)*spin_system.inter.zeeman.matrix{n})*spin_system.inter.basefrqs(n);
                
        end

    end
   
end

% Handle the case when no Zeeman interactions are specified
if (~isfield(inter,'zeeman'))&&(~isfield(inter,'suscept'))
    
    % Warn the user that Zeeman tensors have not been found
    report(spin_system,'WARNING - no Zeeman interactions supplied, magnet frequencies assumed.');
    
    % Fill in magnet frequencies and dipole interaction scaling factors
    for n=1:spin_system.comp.nspins
        spin_system.inter.zeeman.matrix{n}=eye(3)*spin_system.inter.basefrqs(n);
        spin_system.inter.zeeman.ddscal{n}=eye(3);
    end
    
end

% Process giant spin terms
if isfield(inter,'giant')
    
    % Get the array going
    spin_system.inter.giant.coeff=cell(spin_system.comp.nspins,1);
    
    % Loop over spins
    for n=1:spin_system.comp.nspins
        
        % Absorb the coefficients and the Euler angles
        for k=1:numel(inter.giant.coeff{n})
        
            % Get the rotation matrix
            D=wigner(k,inter.giant.euler{n}{k}(1),...
                       inter.giant.euler{n}{k}(2),...
                       inter.giant.euler{n}{k}(3));
               
            % Apply the rotation, convert to angular frequencies and absorb
            spin_system.inter.giant.coeff{n}{k}=2*pi*D*inter.giant.coeff{n}{k}(:);
            
        end
        
    end
   
else
    
    % Warn the user that giant spin terms have not been found
    report(spin_system,'WARNING - no giant spin terms supplied, zeros assumed.');
    
    % Set empty arrays
    spin_system.inter.giant.coeff=cell(spin_system.comp.nspins,1);
   
end

% Absorb chemical subsystems
if isfield(inter,'chem')&&isfield(inter.chem,'parts')
    
    % Assign the data structure
    spin_system.chem.parts=inter.chem.parts;
    
    % Sort spin indices within part specifications
    for n=1:numel(inter.chem.parts)
        inter.chem.parts{n}=sort(inter.chem.parts{n},'ascend');
    end
    
else
    
    % Default is one chemical subsystem
    spin_system.chem.parts={(1:spin_system.comp.nspins)};
    
end

% Absorb concentrations
if isfield(inter,'chem')&&isfield(inter.chem,'concs')
    
    % Assign the data structure
    spin_system.chem.concs=inter.chem.concs;

else
    
    % Default is unit concentration of one subsystem
    if isscalar(spin_system.chem.parts)
        spin_system.chem.concs=1;
    else
        error('concentrations (inter.chem.concs) must be provided.');
    end
    
end

% Absorb exchange reaction rates (linear kinetics)
if isfield(inter,'chem')&&isfield(inter.chem,'rates')
    
    % Assign the data structure
    spin_system.chem.rates=inter.chem.rates;
    
else
    
    % No exchange reactions
    spin_system.chem.rates=[];
    
end

% Absorb magnetization flux rates
if isfield(inter,'chem')&&...
   isfield(inter.chem,'flux_rate')&&...
   isfield(inter.chem,'flux_type')
    
    % Assign the data structure
    spin_system.chem.flux_rate=inter.chem.flux_rate;
    spin_system.chem.flux_type=inter.chem.flux_type;
    
else
    
    % No magnetization fluxes
    spin_system.chem.flux_rate=[];
    spin_system.chem.flux_type='';
    
end

% Report back to the user
summary(spin_system,'chemistry','chemical process summary');

% Order matrices
if isfield(inter,'order_matrix')

    % Absorb the data
    spin_system.inter.order_matrix=inter.order_matrix;
    
    % Inform the user
    report(spin_system,'order matrices supplied by the user:');
    for n=1:numel(spin_system.chem.parts)
        report(spin_system,['  Chemical subsystem ' num2str(n) ':']);
        report(spin_system,['     ' num2str(spin_system.inter.order_matrix{n}(1,:),'%+10.6e  %+10.6e  %+10.6e')]);
        report(spin_system,['     ' num2str(spin_system.inter.order_matrix{n}(2,:),'%+10.6e  %+10.6e  %+10.6e')]);
        report(spin_system,['     ' num2str(spin_system.inter.order_matrix{n}(3,:),'%+10.6e  %+10.6e  %+10.6e')]);
    end
    
else
    
    % Empty cell array
    spin_system.inter.order_matrix={};
    
end
    
% Preallocate coupling tensor array
report(spin_system,'initializing interaction arrays...');
spin_system.inter.coupling.matrix=mat2cell(zeros(3*spin_system.comp.nspins),3*ones(spin_system.comp.nspins,1),...
                                                                            3*ones(spin_system.comp.nspins,1));
% Process coordinates
if isfield(inter,'coordinates')

    % Absorb coordinates into data structure
    spin_system.inter.coordinates=inter.coordinates;
    
    % Report back to the user
    summary(spin_system,'coordinates','atomic coordinates (Angstrom)');
    
    % Process periodic boundary conditions
    if isfield(inter,'pbc')
        
        % Absorb translation vectors into the data structure
        spin_system.inter.pbc=inter.pbc;
        
        % Report back to the user
        summary(spin_system,'pbc','PBC translation vectors (Angstrom)');
        
    else
        
        % Write an empty array
        spin_system.inter.pbc={};
        
        % Report back to the user
        report(spin_system,'periodic boundary conditions not specified, assuming a standalone spin system.');
        
    end
    
    % Call dipolar coupling module
    spin_system=dipolar(spin_system);
    
else
    
    % Warn the user that coordinates have not been found
    report(spin_system,'WARNING - no coordinates given, point dipolar interactions assumed to be zero.');
    
    % Set an empty coordinate array
    spin_system.inter.coordinates=cell(spin_system.comp.nspins,1);
    
    % Set proximity matrix to isolated spins
    spin_system.inter.proxmatrix=speye(spin_system.comp.nspins,spin_system.comp.nspins);
    
end

% Absorb user-specified couplings
if isfield(inter,'coupling')
    
    % Inform the user
    report(spin_system,'processing coupling data...');
    
    % Absorb eigenvalues and Euler angles
    if isfield(inter.coupling,'eigs')
        
        % Get coupling norms
        norms=cellfun(@(x)norm(diag(x),2),inter.coupling.eigs);
        
        % Find significant couplings (Hz at this point)
        [rows,cols,~]=find(norms>spin_system.tols.inter_cutoff);
        
        % Loop over significant couplings
        for n=1:numel(rows)
            
            % Convert Euler angles into DCM
            if isempty(inter.coupling.euler{rows(n),cols(n)})
                S=eye(3,3);
            else
                S=euler2dcm(inter.coupling.euler{rows(n),cols(n)});
            end
            
            % Rotate and convert to rad/s
            component=2*pi*S*diag(inter.coupling.eigs{rows(n),cols(n)})*S';
            
            % Tidy up the numerics
            component=(component+component')/2;
            
            % Add to the total
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+component;
        
        end
        
    end
 
    % Absorb coupling tensors
    if isfield(inter.coupling,'matrix')
        
        % Find significant couplings using 2-norm criterion (Hz at this point)
        [rows,cols,~]=find(cellfun(@(x)norm(x,2),inter.coupling.matrix)>spin_system.tols.inter_cutoff);
        
        % Loop over significant couplings
        for n=1:numel(rows)
            
            % Convert to rad/s and add to the total
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+2*pi*inter.coupling.matrix{rows(n),cols(n)};
        
        end
        
    end
        
    % Absorb scalar couplings
    if isfield(inter.coupling,'scalar')
        
        % Find significant couplings using 2-norm criterion (Hz at this point)
        [rows,cols,~]=find(cellfun(@(x)norm(x,2),inter.coupling.scalar)>spin_system.tols.inter_cutoff);
        
        % Loop over significant couplings
        for n=1:numel(rows)
            
            % Convert to rad/s and add to the total
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+2*pi*eye(3)*inter.coupling.scalar{rows(n),cols(n)};
        
        end
        
    end
    
else
    
    % Warn the user that couplings have not been found
    report(spin_system,'WARNING - no couplings given, zeros assumed.');

end

% Clean up the result
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        
        % Re-check all coupling tensors for significance in 2-norm and relevance
        if (norm(spin_system.inter.coupling.matrix{n,k},2)<2*pi*spin_system.tols.inter_cutoff)||...
           (spin(spin_system.comp.isotopes{n})==0)||(spin(spin_system.comp.isotopes{k})==0)
       
            % Set irrelevant tensors to empty
            spin_system.inter.coupling.matrix{n,k}=[];
            
        end
        
    end
end

% Check for inter-subsystem couplings
for n=1:numel(spin_system.chem.parts)
    for k=1:numel(spin_system.chem.parts)
        if (n~=k)
            coupling_block=spin_system.inter.coupling.matrix(spin_system.chem.parts{n},spin_system.chem.parts{k});
            if ~all(cellfun(@isempty,coupling_block(:)))
                error('couplings detected between spins in different chemical species.');
            end
        end
    end
end

% Report back to the user
summary(spin_system,'couplings','summary of coupling tensors and their anisotropies (Hz)');

% Temperature
if ~isfield(inter,'temperature')
    
    % Default is room temperature
    spin_system.rlx.temperature=298;
    
    % Print a warning to the user
    report(spin_system,'WARNING - spin temperature not specified, assuming 298 Kelvin');
    
else
    
    % Absorb the temperature specified
    spin_system.rlx.temperature=inter.temperature;
    
    % Report back to the user
    report(spin_system,['spin temperature: ' num2str(spin_system.rlx.temperature) ' Kelvin']);
    
end

% Relaxation superoperator
if isfield(inter,'relaxation')
    spin_system.rlx.theories=inter.relaxation;
    for n=1:numel(spin_system.rlx.theories)
        report(spin_system,['relaxation theory ' num2str(n) ': ' spin_system.rlx.theories{n}]);
    end
else
    spin_system.rlx.theories={};
    report(spin_system,'WARNING - no relaxation theory specified');
end

% Rotational correlation times
if isfield(inter,'tau_c')
    spin_system.rlx.tau_c=inter.tau_c;
    for s=1:numel(inter.tau_c)
        report(spin_system,['chemical subsystem ' num2str(s) ', rotational correlation time(s): '...
                                                  num2str(spin_system.rlx.tau_c{s}) ' seconds']);
    end
else
    spin_system.rlx.tau_c={};
end

% Terms to keep in the relaxation superoperator
if ~isempty(spin_system.rlx.theories)
    spin_system.rlx.keep=inter.rlx_keep;
    report(spin_system,['terms to keep in the relaxation superoperator: ' spin_system.rlx.keep]);
end

% The fate of the dynamic frequency shift
if isfield(inter,'rlx_dfs')
    spin_system.rlx.dfs=inter.rlx_dfs;
else
    spin_system.rlx.dfs='ignore';
end
if ~isempty(spin_system.rlx.theories)
    report(spin_system,['action to take on dynamic frequency shifts: ' spin_system.rlx.dfs]);
end

% SRFK correlation time
if isfield(inter,'srfk_tau_c')
    spin_system.rlx.srfk_tau_c=inter.srfk_tau_c;
else
    spin_system.rlx.srfk_tau_c={};
end
if ismember('SRFK',spin_system.rlx.theories)
    for n=1:numel(spin_system.rlx.srfk_tau_c)
        report(spin_system,['SRFK autocorrelation function component ' num2str(n) ':']);
        report(spin_system,['  component weight: ' num2str(spin_system.rlx.srfk_tau_c{n}(1))]);
        report(spin_system,['  correlation time: ' num2str(spin_system.rlx.srfk_tau_c{n}(2)) ' seconds']);
    end
end

% SRFK modulation depths
if isfield(inter,'srfk_mdepth')
    spin_system.rlx.srfk_mdepth=inter.srfk_mdepth;
else
    spin_system.rlx.srfk_mdepth=[];
end
if ismember('SRFK',spin_system.rlx.theories)&&(~isempty(spin_system.rlx.srfk_mdepth))
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            if ~isempty(spin_system.rlx.srfk_mdepth{n,k})
                report(spin_system,['SRFK modulation depth for spins ' ...
                                    num2str(n) ',' num2str(k) ': ' ...
                                    num2str(spin_system.rlx.srfk_mdepth{n,k}) ' Hz']);
            end
        end
    end
end

% SRSK sources
if isfield(inter,'srsk_sources')
    spin_system.rlx.srsk_sources=inter.srsk_sources;
else
    spin_system.rlx.srsk_sources=[];
end

% Equilibrium state
if ~isempty(spin_system.rlx.theories)
    
    % Absorb user setting
    spin_system.rlx.equilibrium=inter.equilibrium;
    
    % Print back the notice
    report(spin_system,['thermalisation method: ' spin_system.rlx.equilibrium]);
    
end

% Global damping
if isfield(inter,'damp_rate')
    spin_system.rlx.damp_rate=inter.damp_rate;
    if isscalar(spin_system.rlx.damp_rate)
        report(spin_system,['global damping rate (Hz): ' num2str(spin_system.rlx.damp_rate,'%+5.3e  ')]);
    else
        report(spin_system,['global damping tensor (Hz): ' num2str(spin_system.rlx.damp_rate(1,:),'%+5.3e  ')]);
        report(spin_system,['                            ' num2str(spin_system.rlx.damp_rate(2,:),'%+5.3e  ')]);
        report(spin_system,['                            ' num2str(spin_system.rlx.damp_rate(3,:),'%+5.3e  ')]);
    end
else
    spin_system.rlx.damp_rate=[];
end

% User-supplied relaxation rates for T1/T2 model
if isfield(inter,'r1_rates')
    spin_system.rlx.r1_rates=inter.r1_rates;
    report(spin_system,'Extended T1/T2 model, summary of R1 relaxation rates:');
    for n=1:numel(spin_system.rlx.r1_rates)
        spin_name=[num2str(n) ' (' spin_system.comp.isotopes{n} ')']; eq_blanks=blanks(numel(spin_name));
        if isscalar(spin_system.rlx.r1_rates{n})
            report(spin_system,['   ' spin_name ' R1 rate (Hz): ' num2str(spin_system.rlx.r1_rates{n},'%+5.3e  ')]);
        else
            report(spin_system,['   ' spin_name ' R1 tensor (Hz): ' num2str(spin_system.rlx.r1_rates{n}(1,:),'%+5.3e  ')]);
            report(spin_system,['   ' eq_blanks '                 ' num2str(spin_system.rlx.r1_rates{n}(2,:),'%+5.3e  ')]);
            report(spin_system,['   ' eq_blanks '                 ' num2str(spin_system.rlx.r1_rates{n}(3,:),'%+5.3e  ')]);
        end
    end
else
    spin_system.rlx.r1_rates={};
end
if isfield(inter,'r2_rates')
    spin_system.rlx.r2_rates=inter.r2_rates;
    report(spin_system,'Extended T1/T2 model, summary of R2 relaxation rates:');
    for n=1:numel(spin_system.rlx.r2_rates)
        spin_name=[num2str(n) ' (' spin_system.comp.isotopes{n} ')']; eq_blanks=blanks(numel(spin_name));
        if isscalar(spin_system.rlx.r2_rates{n})
            report(spin_system,['   ' spin_name ' R2 rate (Hz): ' num2str(spin_system.rlx.r2_rates{n},'%+5.3e  ')]);
        else
            report(spin_system,['   ' spin_name ' R2 tensor (Hz): ' num2str(spin_system.rlx.r2_rates{n}(1,:),'%+5.3e  ')]);
            report(spin_system,['   ' eq_blanks '                 ' num2str(spin_system.rlx.r2_rates{n}(2,:),'%+5.3e  ')]);
            report(spin_system,['   ' eq_blanks '                 ' num2str(spin_system.rlx.r2_rates{n}(3,:),'%+5.3e  ')]);
        end
    end
    report(spin_system,'Extended T1/T2 model: product states will relax at the sum of constituent rates')
else
    spin_system.rlx.r2_rates={};
end

% User-supplied R1 relaxation rates for Lindblad theory
if isfield(inter,'lind_r1_rates')
    spin_system.rlx.lind_r1_rates=inter.lind_r1_rates;
else
    spin_system.rlx.lind_r1_rates=[];
end

% User-supplied R2 relaxation rates for Lindblad theory
if isfield(inter,'lind_r2_rates')
    spin_system.rlx.lind_r2_rates=inter.lind_r2_rates;
else
    spin_system.rlx.lind_r2_rates=[];
end

% User-supplied Weizmann R1e relaxation rates
if isfield(inter,'weiz_r1e')
    spin_system.rlx.weiz_r1e=inter.weiz_r1e;
else
    spin_system.rlx.weiz_r1e=[];
end

% User-supplied Nottingham R1e relaxation rates
if isfield(inter,'nott_r1e')
    spin_system.rlx.nott_r1e=inter.nott_r1e;
else
    spin_system.rlx.nott_r1e=[];
end

% User-supplied Weizmann R1n relaxation rates
if isfield(inter,'weiz_r1n')
    spin_system.rlx.weiz_r1n=inter.weiz_r1n;
else
    spin_system.rlx.weiz_r1n=[];
end

% User-supplied Nottingham R1n relaxation rates
if isfield(inter,'nott_r1n')
    spin_system.rlx.nott_r1n=inter.nott_r1n;
else
    spin_system.rlx.nott_r1n=[];
end

% User-supplied Weizmann R1d relaxation rates
if isfield(inter,'weiz_r1d')
    spin_system.rlx.weiz_r1d=inter.weiz_r1d;
else
    spin_system.rlx.weiz_r1d=[];
end

% User-supplied Weizmann R2e relaxation rates
if isfield(inter,'weiz_r2e')
    spin_system.rlx.weiz_r2e=inter.weiz_r2e;
else
    spin_system.rlx.weiz_r2e=[];
end

% User-supplied Nottingham R2e relaxation rates
if isfield(inter,'nott_r2e')
    spin_system.rlx.nott_r2e=inter.nott_r2e;
else
    spin_system.rlx.nott_r2e=[];
end

% User-supplied Weizmann R2n relaxation rates
if isfield(inter,'weiz_r2n')
    spin_system.rlx.weiz_r2n=inter.weiz_r2n;
else
    spin_system.rlx.weiz_r2n=[];
end

% User-supplied Nottingham R2n relaxation rates
if isfield(inter,'nott_r2n')
    spin_system.rlx.nott_r2n=inter.nott_r2n;
else
    spin_system.rlx.nott_r2n=[];
end

% User-supplied Weizmann R2d relaxation rates
if isfield(inter,'weiz_r2d')
    spin_system.rlx.weiz_r2d=inter.weiz_r2d;
else
    spin_system.rlx.weiz_r2d=[];
end

% Report relaxation rates back to the user
if isfield(inter,'lind_r1_rates')&&isfield(inter,'lind_r2_rates')
    summary(spin_system,'rlx_rates_lindblad','relaxation rates (Hz) for Lindblad theory');
end
if isfield(inter,'nott_r1e')&&isfield(inter,'nott_r2e')&&...
   isfield(inter,'nott_r1n')&&isfield(inter,'nott_r2n')
    summary(spin_system,'rlx_rates_nott','relaxation rates (Hz) for Nottingham DNP theory');
end
if isfield(inter,'weiz_r1e')&&isfield(inter,'weiz_r2e')&&...
   isfield(inter,'weiz_r1n')&&isfield(inter,'weiz_r2n')&&...
   isfield(inter,'weiz_r1d')&&isfield(inter,'weiz_r2d')
    summary(spin_system,'rlx_rates_weiz','relaxation rates (Hz) for Weizmann DNP theory');
end

% Absorb radical recombination parameters
if isfield(inter,'chem')&&isfield(inter.chem,'rp_theory')

    % Absorb theory
    spin_system.chem.rp_theory=inter.chem.rp_theory;
    report(spin_system,['radical recombination theory set to ' spin_system.chem.rp_theory]);

    % Absorb spins
    spin_system.chem.rp_electrons=inter.chem.rp_electrons;
    report(spin_system,['recombining electrons at positions ' num2str(spin_system.chem.rp_electrons)]);

    % Absorb rates
    spin_system.chem.rp_rates=inter.chem.rp_rates;
    report(spin_system,['singlet recombination rate ' num2str(spin_system.chem.rp_rates(1)) ' Hz.']);
    report(spin_system,['triplet recombination rate ' num2str(spin_system.chem.rp_rates(2)) ' Hz.']);

else
    
    spin_system.chem.rp_theory='';
    spin_system.chem.rp_electrons=[];
    spin_system.chem.rp_rates=[];
    
end

% Catch unparsed options
unparsed=fieldnames(sys);
if ~isempty(unparsed)
    for n=1:numel(unparsed)
        report(spin_system,['unrecognised option - ' unparsed{n}]);
    end
    error('unrecognised options in sys');
end

end

% Input validation function
function grumble(sys,inter)

% Check the output switch
if isfield(sys,'output')
    if ~ischar(sys.output)
        error('sys.output must be a character string.');
    end
end

% Check the scratch folder
if isfield(sys,'scratch')
    if ~ischar(sys.scratch)
        error('sys.scratch must be a character string.');
    end
    if ~exist(sys.scratch,'dir')
        error('the specified scratch directory does not exist.');
    end
end

% Check the disable switch
if isfield(sys,'disable')
    if (~iscell(sys.disable))||any(~cellfun(@ischar,sys.disable))
        error('sys.disable must be a cell array of strings.');
    end
    if any(~ismember(sys.disable,{'zte','pt','symmetry','krylov','clean-up','hygiene',...
                                  'dss','expv','trajlevel','merge','colorbar','asyredf'}))
        error('unrecognised switch in sys.disable field.');
    end
end

% Check the enable switch
if isfield(sys,'enable')
    if (~iscell(sys.enable))||any(~cellfun(@ischar,sys.enable))
        error('sys.enable must be a cell array of strings.');
    end
    if any(~ismember(sys.enable,{'gpu','op_cache','xmemlist','greedy','paranoia',...
                                 'cowboy','polyadic','dafuq','prop_cache'}))
        error('unrecognised switch in sys.enable field.');
    end
end

% Check GPU parameters
if isfield(sys,'gpuids')
    
    % Check the specification
    if ~isnumeric(sys.gpuids)||any(mod(sys.gpuids,1)~=0)
        error('sys.gpuids must be a vector of integers.');
    end
    
    % Check if the devices are there
    if any(sys.gpuids>gpuDeviceCount)
        error('GPU device with the specified ID does not exist.');
    end
    
end

% Check isotopes variable
if ~isfield(sys,'isotopes')
    error('sys.isotopes subfield must be present.');
elseif isempty(sys.isotopes)
    error('sys.isotopes cell array cannot be empty.');
elseif ~iscell(sys.isotopes)
    error('sys.isotopes must be a cell array.');
elseif ~all(cellfun(@ischar,sys.isotopes))
    error('all elements of sys.isotopes cell array must be character strings.');
end

% Check labels variable
if isfield(sys,'labels')
    if ~iscell(sys.labels)
        error('sys.labels must be a cell array.');
    elseif isempty(sys.labels)
        error('sys.labels cell array cannot be empty.');
    elseif ~all(cellfun(@ischar,sys.labels))
        error('all elements of sys.labels cell array must be character strings.');
    elseif numel(sys.labels)~=numel(sys.isotopes)
        error('the length of sys.labels must be the same as the length of sys.isotopes.');
    end
end

% Check order matrices
if isfield(inter,'order_matrix')

    % Make sure it's a cell array
    if ~iscell(inter.order_matrix)
        error('inter.order_matrix must be a cell array of 3x3 matrices.');
    end
    
    % Check cell array elements
    for n=1:numel(inter.order_matrix)

        % Make sure they are 3x3 matrices
        if (~isnumeric(inter.order_matrix{n}))||(~ismatrix(inter.order_matrix{n}))||...
           (any(size(inter.order_matrix{n})~=[3 3]))||(~isreal(inter.order_matrix{n}))
            error('elements of inter.order_matrix must be 3x3 matrices of real numbers.');
        end

        % Make sure they are traceless
        if trace(inter.order_matrix{n})>10*eps()
            error('elements of inter.order_matrix must be traceless 3x3 matrices.');
        end

        % Make sure they are symmetric
        if norm(inter.order_matrix{n}-transpose(inter.order_matrix{n}),1)>10*eps()
            error('elements of inter.order_matrix must be symmetric 3x3 matrices.');
        end

    end

    % Make sure it's consistent with chemistry
    if (~isfield(inter,'chem'))||(~isfield(inter.chem,'parts'))
        if numel(inter.order_matrix)~=1
            error('the number of matrices in inter.order_matrix much match the number of chemical subsystems.');
        end
    end
    if isfield(inter,'chem')&&isfield(inter.chem,'parts')&&(numel(inter.order_matrix)~=numel(inter.chem.parts))
        error('the number of matrices in inter.order_matrix much match the number of chemical subsystems.');
    end
    
end

% Check magnet induction
if ~isfield(sys,'magnet')
    error('magnet induction must be specified in sys.magnet varaible.');
else
    if (~isnumeric(sys.magnet))||(~isreal(sys.magnet))||(numel(sys.magnet)~=1)
        error('sys.magnet must be a real number.');
    end
end

% Check Zeeman interactions
if isfield(inter,'zeeman')
    
    % Check eigenvalues / Euler angles specification
    if isfield(inter.zeeman,'eigs')
        
        % Check type
        if (~iscell(inter.zeeman.eigs))||any(~cellfun(@isnumeric,inter.zeeman.eigs))
            error('inter.zeeman.eigs must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Check dimensions
        if numel(inter.zeeman.eigs)~=numel(sys.isotopes)
            error('the number of elements in inter.zeeman.eigs must match the number of spins.');
        end
        
        % Make sure eulers exist
        if ~isfield(inter.zeeman,'euler')
            error('inter.zeeman.euler variable must be set together with inter.zeeman.eigs.');
        end
        
        % Make sure eulers are cells
        if (~iscell(inter.zeeman.euler))||any(~cellfun(@isnumeric,inter.zeeman.euler))
            error('inter.zeeman.euler must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Make sure eulers have the correct length
        if ~all(size(inter.zeeman.eigs)==size(inter.zeeman.euler))
            error('inter.zeeman.eigs and inter.zeeman.euler variables must have the same dimension.');
        end
        
        % Make sure all non-empty elements are real 1x3 vectors
        for n=1:numel(sys.isotopes)
            
            % For eigenvalues
            if ~isempty(inter.zeeman.eigs{n})
                if (~all(size(inter.zeeman.eigs{n})==[1 3]))||...
                   (~isnumeric(inter.zeeman.eigs{n}))||...
                   (~isreal(inter.zeeman.eigs{n}))
                    error('non-empty elements of inter.zeeman.eigs must be real 1x3 vectors.');
                end
            end
            
            % For Euler angles
            if ~isempty(inter.zeeman.euler{n})
                if (~all(size(inter.zeeman.euler{n})==[1 3]))||...
                   (~isnumeric(inter.zeeman.euler{n}))||...
                   (~isreal(inter.zeeman.euler{n}))
                    error('non-empty elements of inter.zeeman.euler must be real 1x3 vectors.');
                end
            end
            
            % For simultaneity
            if (isempty(inter.zeeman.eigs{n})&&(~isempty(inter.zeeman.euler{n})))||...
               (isempty(inter.zeeman.euler{n})&&(~isempty(inter.zeeman.eigs{n})))
                error('inter.zeeman.eigs and inter.zeeman.euler must have identical non-empty cell patterns.');
            end
            
        end
        
    end
    
    % Check matrix specification
    if isfield(inter.zeeman,'matrix')
        
        % Check type
        if ~iscell(inter.zeeman.matrix)
            error('inter.zeeman.matrix must be a cell array of empty matrices or 3x3 matrices.');
        elseif size(inter.zeeman.matrix,1)~=1
            error('inter.zeeman.matrix cell array must have dimension 1 x nspins.');
        elseif ~all(cellfun(@isnumeric,inter.zeeman.matrix))
            error('all elements of inter.zeeman.matrix cell array must be numeric.');
        end 
        
        % Check length
        if numel(inter.zeeman.matrix)~=numel(sys.isotopes)
            error('the number of elements in the inter.zeeman.matrix array should match the number of spins.');
        end
        
        % Make sure all non-empty elements are real 3x3 matrices
        for n=1:numel(sys.isotopes)
            if ~isempty(inter.zeeman.matrix{n})
                if (~all(size(inter.zeeman.matrix{n})==[3 3]))||...
                   (~isnumeric(inter.zeeman.matrix{n}))||...
                   (~isreal(inter.zeeman.matrix{n}))
                    error('non-empty elements of inter.zeeman.matrix must be real 3x3 matrices.');
                end
            end
        end
        
    end
    
    % Check scalars
    if isfield(inter.zeeman,'scalar')
        
        % Check type
        if ~iscell(inter.zeeman.scalar)||any(~cellfun(@isnumeric,inter.zeeman.scalar))
            error('inter.zeeman.scalar must be a cell array of empty matrices or 1x1 matrices.');
        end
        
        % Check length
        if numel(inter.zeeman.scalar)~=numel(sys.isotopes)
            error('the number of elements in the inter.zeeman.scalar array must match the number of spins.');
        end
        
        % Make sure all non-empty elements are real numbers
        for n=1:numel(sys.isotopes)
            if ~isempty(inter.zeeman.scalar{n})
                if (numel(inter.zeeman.scalar{n})~=1)||...
                   (~isnumeric(inter.zeeman.scalar{n}))||...
                   (~isreal(inter.zeeman.scalar{n}))
                    error('non-empty elements of inter.zeeman.scalar must be numbers.');
                end
            end
        end
        
    end
    
end

% Check giant spin parameters
if isfield(inter,'giant')
    
    % Check subfields
    if ~isfield(inter.giant,'coeff')
        error('giant spin model coefficients must be supplied in inter.giant.coeff cell array.');
    end
    if ~isfield(inter.giant,'euler')
        error('giant spin model Euler angles must be supplied in inter.giant.euler cell array.');
    end
    
    % Check the elements
    if numel(inter.giant.coeff)~=numel(sys.isotopes)
        error('the number of elements inter.giant.coeff must be equal to the number of spins.');
    end
    if numel(inter.giant.euler)~=numel(sys.isotopes)
        error('the number of elements inter.giant.euler must be equal to the number of spins.');
    end
    
    % Check the sub-elements
    for n=1:numel(inter.giant.coeff)
        
        % Enforce electrons
        if (~strcmp(sys.isotopes{n}(1),'E'))&&(~isempty(inter.giant.coeff{n}))
            error('giant spin model is only available for electron spins.');
        end
        
        % Enforce Euler angles
        if numel(inter.giant.coeff{n})~=numel(inter.giant.euler{n})
            error(['each spherical rank in the giant spin specification for spin '...
                    num2str(n) ' must have its Euler angles specified.']);
        end
        
        % Enforce sub-element numbers
        for k=1:numel(inter.giant.coeff{n})
            if numel(inter.giant.coeff{n}{k})~=2*k+1
                error(['rank ' num2str(k) ' giant spin ST coefficient specification for spin '...
                        num2str(n) ' must be a vector with ' num2str(2*k+1) ' elements.']);
            end
            if numel(inter.giant.euler{n}{k})~=3
                error(['rank ' num2str(k) ' giant spin ST euler angle specification for spin '...
                        num2str(n) ' must be a vector with three elements.']);
            end
        end
        
    end
    
end

% Check coordinates
if isfield(inter,'coordinates')
    
    % Check type
    if ~iscell(inter.coordinates)||any(~cellfun(@isnumeric,inter.coordinates))
        error('inter.coordinates must be a cell array of empty vectors or 1x3 vectors.');
    end
    
    % Check size
    if numel(inter.coordinates)~=numel(sys.isotopes)
        error('the number of elements in inter.coordinates must match the number of spins.')
    end
    
    % Check contents
    for n=1:numel(sys.isotopes)
        
        % Make sure we have real 3-vectors
        if ~isempty(inter.coordinates{n})
            if (~all(size(inter.coordinates{n})==[1 3]))||...
               (~isnumeric(inter.coordinates{n}))||...
               (~isreal(inter.coordinates{n}))
                error('non-empty elements of inter.coordinates must be real 1x3 vectors.');
            end
        end
        
    end
    
end

% Check periodic boundary conditions
if isfield(inter,'pbc')
    
    % Check type
    if ~iscell(inter.pbc)||any(~cellfun(@isvector,inter.pbc))
        error('inter.pbc must be a cell array of row vectors.');
    end
    
    % Check element numbers
    if ~ismember(numel(inter.pbc),[0 1 2 3])
        error('inter.pbc cell array must have zero, one, two or three row vectors as elements.');
    end
    
    % Check vector dimensions
    for n=1:numel(inter.pbc)
        if (~all(size(inter.pbc{n})==[1 3]))||(~isreal(inter.pbc{n}))
            error('all elements of inter.pbc.cell array must be row vectors with three real elements.');
        end
    end
    
    % Check linear independence
    if ((numel(inter.pbc)==2)&&rank([inter.pbc{1}; inter.pbc{2}],1e-3)<2)||...
       ((numel(inter.pbc)==3)&&rank([inter.pbc{1}; inter.pbc{2}; inter.pbc{3}],1e-3)<3)
        error('the vectors supplied in inter.pbc must be lineary independent.');
    end
    
end
            
% Check couplings
if isfield(inter,'coupling')
    
    % Check eigenvalues / Euler angles specification
    if isfield(inter.coupling,'eigs')
        
        % Check type
        if (~iscell(inter.coupling.eigs))||any(any(~cellfun(@isnumeric,inter.coupling.eigs)))
            error('inter.coupling.eigs must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.eigs)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.eigs must match the number of spins.');
        end
        
        % Check quadratic couplings
        for n=1:numel(sys.isotopes)
            [~,mult]=spin(sys.isotopes{n});
            if (norm(inter.coupling.eigs{n,n},1)>0)&&(mult<3)
                error('quadratic couplings cannot be specified for spin-1/2 particles.');
            elseif abs(sum(inter.coupling.eigs{n,n}))>1e-6
                error('quadratic couplings cannot have a non-zero trace.');
            end
        end
        
        % Make sure eulers exist
        if ~isfield(inter.coupling,'euler')
            error('inter.coupling.euler array must be set together with inter.coupling.eigs.');
        end
        
        % Make sure eulers are cells
        if ~iscell(inter.coupling.euler)||any(any(~cellfun(@isnumeric,inter.coupling.euler)))
            error('inter.coupling.euler must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Make sure eulers have the correct length
        if ~all(size(inter.coupling.eigs)==size(inter.coupling.euler))
            error('inter.coupling.eigs and inter.coupling.euler arrays must have the same dimension.');
        end
        
        % Make sure all non-empty elements are real 1x3 vectors
        for n=1:numel(sys.isotopes)
            for k=1:numel(sys.isotopes)
            
                % For eigenvalues
                if ~isempty(inter.coupling.eigs{n,k})
                    if (~all(size(inter.coupling.eigs{n,k})==[1 3]))||...
                       (~isnumeric(inter.coupling.eigs{n,k}))||...
                       (~isreal(inter.coupling.eigs{n,k}))
                        error('non-empty elements of inter.coupling.eigs must be real 1x3 vectors.');
                    end
                end
                
                % For Euler angles
                if ~isempty(inter.coupling.euler{n,k})
                    if (~all(size(inter.coupling.euler{n,k})==[1 3]))||...
                       (~isnumeric(inter.coupling.euler{n,k}))||...
                       (~isreal(inter.coupling.euler{n,k}))
                        error('non-empty elements of inter.coupling.euler must be real 1x3 vectors.');
                    end
                end
                
                % For simultaneity
                if (isempty(inter.coupling.eigs{n,k})&&(~isempty(inter.coupling.euler{n,k})))||...
                   (isempty(inter.coupling.euler{n,k})&&(~isempty(inter.coupling.eigs{n,k})))
                    error('inter.coupling.eigs and inter.coupling.euler must have identical non-empty cell patterns.');
                end
                
            end
        end
        
    end
    
    % Check matrix specification
    if isfield(inter.coupling,'matrix')
        
        % Check type
        if ~iscell(inter.coupling.matrix)||any(any(~cellfun(@isnumeric,inter.coupling.matrix)))
            error('inter.coupling.matrix must be a cell array of empty matrices or 3x3 matrices.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.matrix)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.matrix cell array should match the number of spins.');
        end
        
        % Check quadratic couplings
        for n=1:numel(sys.isotopes)
            [~,mult]=spin(sys.isotopes{n});
            if (norm(inter.coupling.matrix{n,n},1)>0)&&(mult<3)
                error('quadratic couplings cannot be specified for spin-1/2 particles.');
            elseif abs(trace(inter.coupling.matrix{n,n}))>1e-6
                error('quadratic couplings cannot have a non-zero trace.');
            end
        end
        
        % Make sure all non-empty elements are real 3x3 matrices
        for n=1:numel(sys.isotopes)
            for k=1:numel(sys.isotopes)
                if ~isempty(inter.coupling.matrix{n,k})
                    if (~all(size(inter.coupling.matrix{n,k})==[3 3]))||...
                       (~isnumeric(inter.coupling.matrix{n,k}))||...
                       (~isreal(inter.coupling.matrix{n,k}))
                        error('non-empty elements of inter.coupling.matrix must be real 3x3 matrices.');
                    end
                end
            end
        end
        
    end
    
    % Check scalars
    if isfield(inter.coupling,'scalar')
        
        % Check type
        if ~iscell(inter.coupling.scalar)||any(any(~cellfun(@isnumeric,inter.coupling.scalar)))
            error('inter.coupling.scalar must be a cell array of empty matrices or 1x1 matrices.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.scalar)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.scalar array must match the number of spins.');
        end
        
        % Make sure all non-empty elements are real numbers
        if any(nonzeros(cellfun(@numel,inter.coupling.scalar))~=1)||...
           any(nonzeros(~cellfun(@isnumeric,inter.coupling.scalar)))||...
           any(nonzeros(~cellfun(@isreal,inter.coupling.scalar)))
            error('non-empty elements of inter.coupling.scalar must be real numbers.');
        end
        
        % Disallow quadratic scalar couplings
        for n=1:numel(sys.isotopes)
            if norm(inter.coupling.scalar{n,n},1)~=0
                error('scalar couplings cannot be quadratic.');
            end
        end
         
    end
    
end

% Check temperature - negative and complex values allowed
if isfield(inter,'temperature')&&(~isempty(inter.temperature))
    
    % Check type and dimension
    if (~isnumeric(inter.temperature))||(numel(inter.temperature)~=1)
        error('inter.tempearture must be a number.');
    end
    
end

% Check relaxation theory
if isfield(inter,'relaxation')
    
    % Check type
    if (~iscell(inter.relaxation))||any(~cellfun(@ischar,inter.relaxation))
        error('inter.relaxation must be a cell array of strings.');
    end
    
    % Check relaxation theories
    if ~all(ismember(inter.relaxation,{'damp','t1_t2','redfield','lindblad',...
                                       'nottingham','weizmann','SRFK','SRSK'}))
        error('unrecognised relaxation theory specification.');
    end
    
    % Enforce term retention policy specification
    if ~isfield(inter,'rlx_keep')
        error('relaxation superoperator term retention policy must be specified in inter.rlx_keep field.');
    end
    
    % Enforce relaxation destination
    if ~isfield(inter,'equilibrium')
        error('relaxation destination must be specified in inter.equilibrium variable.');
    end
    
    % Enforce correlation time with Redfield theory
    if ismember('redfield',inter.relaxation)&&(~isfield(inter,'tau_c'))
        error('correlation time(s) must be specified with Redfield theory.');
    end
    
    % Enforce damping rate with isotropic damping
    if ismember('damp',inter.relaxation)&&(~isfield(inter,'damp_rate'))
        error('damping rate or tensor must be specified with non-selective damping.');
    end
  
    % Enforce R1 and R2 rates with T1,T2 approximation
    if ismember('t1_t2',inter.relaxation)&&((~isfield(inter,'r1_rates'))||(~isfield(inter,'r2_rates')))
        error('R1 and R2 rates must be specified with extended T1,T2 relaxation theory.');
    end
    
    % Enforce R1 and R2 rates with Lindblad theory
    if ismember('lindblad',inter.relaxation)&&((~isfield(inter,'lind_r1_rates'))||(~isfield(inter,'lind_r2_rates')))
        error('R1 and R2 rates must be specified with Lindblad relaxation theory.');
    end
    
    % Enforce R1e, R2e, R1n and R2n rates with Nottingham DNP theory
    if ismember('nottingham',inter.relaxation)&&...
       ((~isfield(inter,'nott_r1e'))||(~isfield(inter,'nott_r2e'))||...
        (~isfield(inter,'nott_r1n'))||(~isfield(inter,'nott_r2n')))
        error(['R1e, R2e, R1n and R2n rates must be specified with '...
               'Nottingham DNP relaxation theory.']);
    end
    
    % Enforce R1e, R2e, R1n and R2n rates with Weizmann DNP theory
    if ismember('weizmann',inter.relaxation)&&...
       ((~isfield(inter,'weiz_r1e'))||(~isfield(inter,'weiz_r2e'))||...
        (~isfield(inter,'weiz_r1n'))||(~isfield(inter,'weiz_r2n')))
        error(['R1e, R2e, R1n and R2n rates must be specified with '...
               'Weizmann DNP relaxation theory.']);
    end
    
    % Enforce R1d and R2d with Weizmann DNP theory
    if ismember('weizmann',inter.relaxation)&&...
       ((~isfield(inter,'weiz_r1d'))||(~isfield(inter,'weiz_r2d')))
        error('R1d and R2d rates must be specified with Weizmann DNP relaxation theory.');
    end
    
    % Enforce two electrons with Nottingham DNP theory
    if ismember('nottingham',inter.relaxation)&&(nnz(strcmp('E',sys.isotopes))~=2)
        error('Nottingham DNP relaxation theory requires two electrons.');
    end
    
    % Enforce correlation time with SRFK
    if ismember('SRFK',inter.relaxation)&&(~isfield(inter,'srfk_tau_c'))
        error('SRFK requires modulation correlation time to be specified in inter.srfk_tau_c variable.');
    end
    
    % Enforce modulation depth with SRFK
    if ismember('SRFK',inter.relaxation)&&(~isfield(inter,'srfk_mdepth'))
        error('SRFK requires modulation depths to be specified in inter.srfk_mdepth variable.');
    end

    % Insist on source list with SRSK
    if ismember('SRSK',inter.relaxation)&&(~isfield(inter,'srsk_sources'))
        error('SRSK requires source spin list in inter.srsk_sources variable.');
    end
    
end

% Check rotational correlation time
if isfield(inter,'tau_c')
    
    if ~iscell(inter.tau_c)
        error('inter.tau_c must be a cell array with the number of elements matching the number of chemical species.');
    end
    
    % Inspect the elements
    for s=1:numel(inter.tau_c)
        
        % Check type and dimension
        if (~isnumeric(inter.tau_c{s}))||(numel(inter.tau_c{s})>3)||(numel(inter.tau_c{s})==0)
            error('elements of inter.tau_c must be numeric vectors of size 1, 2 or 3.');
        end
        
        % Check value
        if (~isreal(inter.tau_c{s}))||any(inter.tau_c{s}<0)
            error('elements of inter.tau_c must be non-negative and real.');
        end
        
        % Enforce Redfield theory if tau_c is specified
        if (~isfield(inter,'relaxation'))||(~ismember('redfield',inter.relaxation))
            error('inter.tau_c requires Redfield relaxation theory.');
        end
        
    end
    
    % Match the chemistry
    if isfield(inter,'chem')&&isfield(inter.chem,'parts')
        if numel(inter.chem.parts)~=numel(inter.tau_c)
            error('the number of elements in inter.tau_c must match the number of chemical species.');
        end
    elseif numel(inter.tau_c)~=1
        error('a system with a single chemical species can only have one element in inter.tau_c cell array.');
    end
    
end

% Check SRFK correlation time
if isfield(inter,'srfk_tau_c')
    
    % Check type and dimension
    if ~iscell(inter.srfk_tau_c)
        error('inter.srfk_tau_c must be a cell array of 2-element vectors.');
    end
    for n=1:numel(inter.srfk_tau_c)
        if (~isnumeric(inter.srfk_tau_c{n}))||...
           (~isreal(inter.srfk_tau_c{n}))||...
           (numel(inter.srfk_tau_c{n})~=2)
            error('inter.srfk_tau_c must be a cell array of 2-element vectors.');
        end
    end
    
    % Enforce SRFK theory if srfk_tau_c is specified
    if (~isfield(inter,'relaxation'))||(~ismember('SRFK',inter.relaxation))
        error('inter.srfk_tau_c requires SRFK relaxation theory.');
    end
    
end

% Check term retention
if isfield(inter,'rlx_keep')
    
    % Check type
    if ~ischar(inter.rlx_keep)
        error('inter.rlx_keep must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.rlx_keep,{'diagonal','kite','secular','labframe'})
        error('allowed values for inter.rlx_keep are ''diagonal'', ''kite'', ''secular'' and ''labframe''.');
    end
    
end

% Check dynamic frequency shift retention
if isfield(inter,'rlx_dfs')
    
    % Check type
    if ~ischar(inter.rlx_dfs)
        error('inter.rlx_dfs must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.rlx_dfs,{'keep','ignore'})
        error('allowed values for inter.rlx_dfs are ''keep'' and ''ignore''.');
    end
    
end
    
% Check SRFK modulation depths
if isfield(inter,'srfk_mdepth')

    % Check type
    if (~iscell(inter.srfk_mdepth))||(~all(cellfun(@isnumeric,inter.srfk_mdepth(:))))
        error('inter.srfk_mdepth must be a cell array of empty matrices or scalars.');
    end
    
    % Check dimensions
    if ~all(size(inter.srfk_mdepth)==[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.srfk_mdepth array should match the number of spins.');
    end
    
    % Make sure all non-empty elements are non-negative real numbers
    [rows,cols]=find(~cellfun(@isempty,inter.srfk_mdepth));
    for n=1:numel(rows)
        if (~isreal(inter.srfk_mdepth{rows(n),cols(n)}))||(inter.srfk_mdepth{rows(n),cols(n)}<0)
            error('non-empty elements of inter.srfk_mdepth must be non-negative real numbers.');
        end
    end
    
    % Disallow quadratic scalar couplings
    for n=1:numel(sys.isotopes)
        if norm(inter.srfk_mdepth{n,n},1)~=0
            error('scalar couplings cannot be quadratic.');
        end
    end
    
    % Enforce SRFK if SRFK modulation depths are specified
    if (~isfield(inter,'relaxation'))||(~ismember('SRFK',inter.relaxation))
        error('inter.srfk_mdepth requires SRFK relaxation theory.');
    end
    
end

% Check equilibrium switch
if isfield(inter,'equilibrium')
    
    % Check type
    if ~ischar(inter.equilibrium)
        error('inter.equilibrium must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.equilibrium,{'zero','IME','dibari'})
        error('allowed values for inter.equilibrium are ''zero'', ''IME'' and ''dibari''.');
    end
    
    % Enforce temperature
    if ismember(inter.equilibrium,{'IME','dibari'})&&...
       (~isfield(inter,'temperature'))
        error('inter.temperature is required when relaxation has a target.');
    end
    
end

% Check damping rate
if isfield(inter,'damp_rate')

    % Enforce isotropic damping if damp_rate is specified
    if ~ismember('damp',inter.relaxation)
        error('inter.damp_rate can only be specified with damp relaxation theory.');
    end

    % Check type
    if (~isnumeric(inter.damp_rate))||(~isreal(inter.damp_rate))||...
       ((~isscalar(inter.damp_rate))&&(~all(size(inter.damp_rate)==3)))
        error('inter.damp_rate must be a real scalar or a 3x3 matrix.');
    end
    
    % Check dimensions
    if ~all(eig(inter.damp_rate,'vector')>=0)
        error('inter.damp_rate must have non-negative eigenvalues.');
    end
    
end

% Check R1 rates for T1/T2
if isfield(inter,'r1_rates')
    
    % Check type and count
    if ~iscell(inter.r1_rates)
        error('inter.r1_rates must be a cell array of scalars or 3x3 matrices.')
    end
    if numel(inter.r1_rates)~=numel(sys.isotopes)
        error('element count in inter.r1_rates must match the number of spins.');
    end

    % Check elements
    for n=1:numel(inter.r1_rates)
        if (~isnumeric(inter.r1_rates{n}))||(~isreal(inter.r1_rates{n}))||...
           ((~isscalar(inter.r1_rates{n}))&&(~all(size(inter.r1_rates{n})==[3 3])))
            error('elements of inter.r1_rates must be real scalars or 3x3 matrices.');
        end
        if norm(inter.r1_rates{n}-inter.r1_rates{n}',2)>1e-6*norm(inter.r1_rates{n},2)
            error('3x3 matrices in inter.r1_rates must be symmetric.');
        end
        if any(eig(inter.r1_rates{n},'vector')<0)
            error('elements of inter.r1_rates must be non-negative definite.');
        end
    end
    
    % Enforce T1,T2 theory if inter.r1_rates rates are specified
    if ~ismember('t1_t2',inter.relaxation)
        error('inter.r1_rates can only be specified with T1,T2 relaxation theory.');
    end
    
end

% Check R2 rates for T1/T2
if isfield(inter,'r2_rates')
    
    % Check type and count
    if ~iscell(inter.r2_rates)
        error('inter.r2_rates must be a cell array of scalars or 3x3 matrices.')
    end
    if numel(inter.r2_rates)~=numel(sys.isotopes)
        error('element count in inter.r2_rates must match the number of spins.');
    end

    % Check elements
    for n=1:numel(inter.r2_rates)
        if (~isnumeric(inter.r2_rates{n}))||(~isreal(inter.r2_rates{n}))||...
           ((~isscalar(inter.r2_rates{n}))&&(~all(size(inter.r2_rates{n})==[3 3])))
            error('elements of inter.r2_rates must be real scalars or 3x3 matrices.');
        end
        if norm(inter.r2_rates{n}-inter.r2_rates{n}',2)>1e-6*norm(inter.r2_rates{n},2)
            error('3x3 matrices in inter.r2_rates must be symmetric.');
        end
        if any(eig(inter.r2_rates{n},'vector')<0)
            error('elements of inter.r2_rates must be non-negative definite.');
        end
    end
    
    % Enforce T1,T2 theory if inter.r2_rates rates are specified
    if ~ismember('t1_t2',inter.relaxation)
        error('inter.r2_rates can only be specified with T1,T2 relaxation theory.');
    end
    
end

% Check R1 rates for Lindblad
if isfield(inter,'lind_r1_rates')
    
    % Check type and values
    if (~isvector(inter.lind_r1_rates))||(~isnumeric(inter.lind_r1_rates))||...
         any(inter.lind_r1_rates<0)||any(~isreal(inter.lind_r1_rates))
        error('inter.lind_r1_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.lind_r1_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.lind_r1_rates must be equal to the number of spins.');
    end
    
    % Enforce Lindblad theory if inter.lind_r1_rates rates are specified
    if ~ismember('lindblad',inter.relaxation)
        error('inter.lind_r1_rates can only be specified with Lindblad relaxation theory.');
    end
    
end

% Check R2 rates for Lindblad
if isfield(inter,'lind_r2_rates')
    
    % Check type and values
    if (~isvector(inter.lind_r2_rates))||(~isnumeric(inter.lind_r2_rates))||...
         any(inter.lind_r2_rates<0)||any(~isreal(inter.lind_r2_rates))
        error('inter.lind_r2_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.lind_r2_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.lind_r2_rates must be equal to the number of spins.');
    end
    
    % Enforce Lindblad theory if inter.lind_r2_rates rates are specified
    if ~ismember('lindblad',inter.relaxation)
        error('inter.lind_r2_rates can only be specified with Lindblad relaxation theory.');
    end
    
end

% Check R1e rate for Weizmann DNP theory
if isfield(inter,'weiz_r1e')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r1e))||any(inter.weiz_r1e<0)||any(~isreal(inter.weiz_r1e))
        error('inter.weiz_r1e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r1e)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1e is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r1e can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R1e rate for Nottingham DNP theory
if isfield(inter,'nott_r1e')
    
    % Check type and value
    if (~isnumeric(inter.nott_r1e))||any(inter.nott_r1e<0)||any(~isreal(inter.nott_r1e))
        error('inter.nott_r1e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r1e)~=1
        error('inter.nott_r1e must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r1e can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R2e rate for Weizmann DNP theory
if isfield(inter,'weiz_r2e')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r2e))||any(inter.weiz_r2e<0)||any(~isreal(inter.weiz_r2e))
        error('inter.weiz_r2e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r2e)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r2e is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r2e can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2e rate for Nottingham DNP theory
if isfield(inter,'nott_r2e')
    
    % Check type and value
    if (~isnumeric(inter.nott_r2e))||any(inter.nott_r2e<0)||any(~isreal(inter.nott_r2e))
        error('inter.nott_r2e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r2e)~=1
        error('inter.nott_r2e must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r2e can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R1n rate for Weizmann DNP theory
if isfield(inter,'weiz_r1n')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r1n))||any(inter.weiz_r1n<0)||any(~isreal(inter.weiz_r1n))
        error('inter.weiz_r1n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r1n)~=1
        error('inter.weiz_r1n must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1n is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r1n can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R1n rate for Nottingham DNP theory
if isfield(inter,'nott_r1n')
    
    % Check type and value
    if (~isnumeric(inter.nott_r1n))||any(inter.nott_r1n<0)||any(~isreal(inter.nott_r1n))
        error('inter.nott_r1n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r1n)~=1
        error('inter.nott_r1n must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1n is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r1n can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R2n rate for Weizmann DNP theory
if isfield(inter,'weiz_r2n')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r2n))||any(inter.weiz_r2n<0)||any(~isreal(inter.weiz_r2n))
        error('inter.weiz_r2n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r2n)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r2n is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r2n can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2n rate for Nottingham DNP theory
if isfield(inter,'nott_r2n')
    
    % Check type and value
    if (~isnumeric(inter.nott_r2n))||any(inter.nott_r2n<0)||any(~isreal(inter.nott_r2n))
        error('inter.nott_r2n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r2n)~=1
        error('inter.nott_r2n must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r2n can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R1d rates for Weizmann DNP theory
if isfield(inter,'weiz_r1d')
    
    % Check type and values
    if (~isnumeric(inter.weiz_r1d))||any(inter.weiz_r1d(:)<0)||(~isreal(inter.weiz_r1d))
        error('inter.weiz_r1d must be a matrix of non-negative real numbers.');
    end
    
    % Check dimension
    if any(size(inter.weiz_r1d)~=[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.weiz_r1d must be equal to the number of spins.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1d rates are specified
    if ~strcmp(inter.relaxation,'weizmann')
        error('inter.weiz_r1d can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2d rates
if isfield(inter,'weiz_r2d')
    
    % Check type and values
    if (~isnumeric(inter.weiz_r2d))||any(inter.weiz_r2d(:)<0)||(~isreal(inter.weiz_r2d))
        error('inter.weiz_r2d must be a matrix of non-negative real numbers.');
    end
    
    % Check dimension
    if any(size(inter.weiz_r2d)~=[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.weiz_r2d must be equal to the number of spins.');
    end
    
    % Enforce Weizmann theory if R2d rates are specified
    if ~strcmp(inter.relaxation,'weizmann')
        error('inter.weiz_r2d can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check chemical kinetics
if isfield(inter,'chem')
    
    % If rates are provided, insist on parts and concentrations
    if isfield(inter.chem,'rates')&&(~isfield(inter.chem,'parts'))
        error('subsystem identifiers (inter.chem.parts) must be provided.');
    elseif isfield(inter.chem,'rates')&&(~isfield(inter.chem,'concs'))
        error('initial concentrations (inter.chem.concs) must be provided.');
    end
    
    % Check chemical species specification
    if isfield(inter.chem,'parts')

        % Basic type checks
        if ~iscell(inter.chem.parts)||(~all(cellfun(@isvector,inter.chem.parts)))
            error('inter.chem.parts must be a cell array of vectors.');
        end

        % Chemiscal subsystem specification
        for n=1:numel(inter.chem.parts)

            % Subsystem partitioning check
            for k=1:numel(inter.chem.parts)
                if (n~=k)&&(~isempty(intersect(inter.chem.parts{n},inter.chem.parts{k})))
                    error('a given spin can only belong to one chemical subsystem in inter.chem.parts variable.');
                end
            end

            % Spin indexing check
            if any(inter.chem.parts{n}<1)||any(inter.chem.parts{n}>numel(sys.isotopes))||...
               any(mod(inter.chem.parts{n},1)~=0)||(numel(unique(inter.chem.parts{n}))~=numel(inter.chem.parts{n}))
                error('elements of inter.chem.parts must be vectors of unique positive integers not exceeding the total number of spins.');
            end

            % If rates are given, insist on first order exchange
            if isfield(inter.chem,'rates')
                for k=1:numel(inter.chem.parts)
                    if numel(inter.chem.parts{n})~=numel(inter.chem.parts{k})
                        error('exchange mode: all chemical subsystems must have the same number of spins.');
                    end
                end
                for k=1:numel(inter.chem.parts)
                    for m=1:numel(inter.chem.parts{n})
                        if ~strcmp(sys.isotopes{inter.chem.parts{n}(m)},sys.isotopes{inter.chem.parts{k}(m)})
                            error('exchange mode: isotope sequences in all chemical subsystems must be the same.');
                        end
                    end
                end
            end

        end

    end
    
    % Check reaction rate matrix
    if isfield(inter.chem,'rates')
        if (~isnumeric(inter.chem.rates))||(~isreal(inter.chem.rates))||...
           (size(inter.chem.rates,1)~=size(inter.chem.rates,2))
            error('inter.chem.rates must be a real square matrix.');
        end
        if any(size(inter.chem.rates)~=numel(inter.chem.parts))
            error('both dimensions of inter.chem.rates matrix must be equal to the number of chemical subsystems.');
        end
        if ~all(abs(sum(inter.chem.rates,1))<10*eps('double'))
            error('inter.chem.rates violates conservation of matter: column sums must be zero.');
        end
    end
    
    % Check initial concentrations
    if isfield(inter.chem,'concs')
        if (~isnumeric(inter.chem.concs))||(~isreal(inter.chem.concs))||any(inter.chem.concs(:)<0)
            error('inter.chem.concs must be a vector of non-negative real numbers.');
        end
        if numel(inter.chem.concs)~=numel(inter.chem.parts)
            error('the number of initial concentrations must be equal to the number of chemical species.');
        end
    end
    
    % Check flux specifications
    if isfield(inter.chem,'flux_rate')&&(~isfield(inter.chem,'flux_type'))
        error('flux type (inter.chem.flux_type) must be provided.');
    elseif isfield(inter.chem,'flux_type')&&(~isfield(inter.chem,'flux_rate'))
        error('flux rates (inter.chem.flux_rate) must be provided.');
    end
    
    % Check flux rate matrix
    if isfield(inter.chem,'flux_rate')
        if (~isnumeric(inter.chem.flux_rate))||(~isreal(inter.chem.flux_rate))||...
           (size(inter.chem.flux_rate,1)~=size(inter.chem.flux_rate,2))
            error('inter.chem.flux_rate must be a real square matrix.');
        end
        if any(size(inter.chem.flux_rate)~=numel(sys.isotopes))
            error('both dimensions of inter.chem.flux_rate matrix must be equal to the number of spins.');
        end
    end
    
    % Check flux type
    if isfield(inter.chem,'flux_type')
        if ~ischar(inter.chem.flux_type)
            error('inter.chem.flux_type must be a character string.');
        end
        if ~ismember(inter.chem.flux_type,{'intermolecular','intramolecular'})
            error('incorrect flux type specification.');
        end
    end
    
    % Check radical pair kinetics
    if isfield(inter.chem,'rp_theory')
        if ~ischar(inter.chem.rp_theory)
            error('inter.chem.rp_theory must be a string.');
        end
        if ~ismember(inter.chem.rp_theory,{'haberkorn','jones-hore','exponential'})
            error('allowed values for inter.chem.rp_theory are ''exponential'', ''haberkorn'' and ''jones-hore''.');
        end
        if (~isfield(inter.chem,'rp_electrons'))||(~isfield(inter.chem,'rp_rates'))
            error('inter.chem.rp_electrons and inter.chem.rp_rates must be specified alongside inter.chem.rp_theory parameter.');
        end
    end
    if isfield(inter.chem,'rp_electrons')
        if (~isfield(inter.chem,'rp_theory'))||(~isfield(inter.chem,'rp_rates'))
            error('inter.chem.rp_theory and inter.chem.rp_rates must be specified alongside inter.chem.rp_electrons parameter.');
        end
        if (~isnumeric(inter.chem.rp_electrons))||(numel(inter.chem.rp_electrons)~=2)||...
           any(mod(inter.chem.rp_electrons,1)~=0)||any(inter.chem.rp_electrons<1)
            error('inter.chem.rp_electrons must be a vector of two positive integers.');
        end
        if any(inter.chem.rp_electrons>numel(sys.isotopes))||any(~cellfun(@(x)strcmp(x(1),'E'),sys.isotopes(inter.chem.rp_electrons)))
            error('at least one of the elements of inter.chem.rp_electrons does not refer to an electron.');
        end
    end
    if isfield(inter.chem,'rp_rates')
        if (~isfield(inter.chem,'rp_theory'))||(~isfield(inter.chem,'rp_electrons'))
            error('inter.chem.rp_theory and inter.chem.rp_electrons must be specified alongside inter.chem.rp_rates parameter.');
        end
        if (~isnumeric(inter.chem.rp_rates))||(numel(inter.chem.rp_rates)~=2)||(~isreal(inter.chem.rp_rates))||any(inter.chem.rp_rates<0)
            error('inter.chem.rp_rates must be a vector of two non-negative real numbers.');
        end
    end
    
end

end

% Those who beat their swords into plowshares will till 
% the soil for those who did not.
%
% Benjamin Franklin

