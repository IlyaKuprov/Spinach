% Comparison of powder-averaged sliced and Fokker-Planck MAS
% evolution for a 13C spin with CSA and transverse RF. Both routes
% use the same weighted Lebedev crystallite grid.
%
% Syntax: mas_fplanck_powder()
%
% Checks rotor rank, slice count, rotor phase, and powder grid.
%
% Calculation time: minutes
%
% talos@spindynamics.org

function mas_fplanck_powder()

% System specification
sys.magnet=9.4;
sys.isotopes={'13C'};
sys.parallel={'processes',1};
sys.parprops={};
inter.zeeman.eigs={[-120 -25 145]};
inter.zeeman.euler={[0.4 0.7 0.2]};

% Basis set
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spin system
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system.sys.output='hush';

% MAS and RF parameters
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='leb_2ang_rank_5';
parameters.spins={'13C'};
parameters.offset=0;
parameters.rframes={{'13C',1}};
parameters.rho0=state(spin_system,'L+','13C');
parameters.coil=parameters.rho0;
parameters.ref_norm=parameters.coil'*parameters.rho0;
parameters.rf_op=operator(spin_system,'Lx','13C');
parameters.rf_amp=2*pi*3500;
parameters.duration=1/parameters.rate;
parameters.serial=true;
parameters.verbose=0;

% Fokker-Planck rotor-rank convergence
rotor_ranks=[6 8];
fp_sig=zeros(size(rotor_ranks));
for rank_idx=1:numel(rotor_ranks)
    parameters.max_rank=rotor_ranks(rank_idx);
    fp_sig(rank_idx)=singlerot(spin_system,@fp_powder_signal,...
                               parameters,'labframe');
    fprintf('Powder FP rank %d: %.9g%+.9gi\n',rotor_ranks(rank_idx),...
            real(fp_sig(rank_idx)),imag(fp_sig(rank_idx)));
end

% Sliced propagation and rotor-phase convergence
slice_counts=[33 65];
phase_counts=[7 13];
sl_sig=zeros(numel(slice_counts),numel(phase_counts));
for slice_idx=1:numel(slice_counts)
    for phase_idx=1:numel(phase_counts)
        sl_sig(slice_idx,phase_idx)=sliced_powder_signal(...
            spin_system,parameters,slice_counts(slice_idx),...
            phase_counts(phase_idx));
        fprintf('Powder slices %d phases %d: %.9g%+.9gi\n',...
                slice_counts(slice_idx),phase_counts(phase_idx),...
                real(sl_sig(slice_idx,phase_idx)),...
                imag(sl_sig(slice_idx,phase_idx)));
    end
end

% Powder quadrature convergence
parameters.grid='leb_2ang_rank_11';
fp_fine=singlerot(spin_system,@fp_powder_signal,...
                  parameters,'labframe');
grid_step=abs(fp_fine-fp_sig(end));

% Comparison of the two simulation routes
fp_step=abs(fp_sig(end)-fp_sig(end-1));
sl_step=abs(sl_sig(end,end)-sl_sig(end-1,end));
phase_step=abs(sl_sig(end,end)-sl_sig(end,end-1));
route_gap=abs(fp_sig(end)-sl_sig(end,end));
fprintf('Powder CSA: FP step %.6g, slice step %.6g, ',fp_step,sl_step);
fprintf('phase step %.6g, grid step %.6g, gap %.6g\n',...
        phase_step,grid_step,route_gap);
tolerance=0.002;
assert(all(isfinite([fp_step sl_step phase_step grid_step ...
                     route_gap]))&&...
       fp_step<tolerance/2&&sl_step<tolerance/2&&...
       phase_step<tolerance/2&&grid_step<tolerance&&route_gap<tolerance,...
       'Powder-averaged MAS routes did not converge.');
fprintf('MAS_FPLANCK_POWDER_SUCCESS\n');

end

% Build the sliced propagator for each crystallite and rotor start phase
function signal=sliced_powder_signal(spin_system,parameters,...
                                     slice_count,phase_count)

% Load the common powder grid and configure single-crystal rotor stacks
powder_grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',...
                   parameters.grid),'alphas','betas','gammas','weights');
parameters.grid='single_crystal';
parameters.masframe='rotor';
parameters.max_rank=(slice_count-1)/2;

% Average rotor start phases for each weighted crystallite
signal=0;
for crystal_idx=1:numel(powder_grid.weights)
    crystal_signal=0;
    for phase_idx=1:phase_count
        rotor_phase=2*pi*(phase_idx-1)/phase_count;

        % Build the phase-shifted midpoint rotor stack
        parameters.orientation=[powder_grid.alphas(crystal_idx)+...
                                rotor_phase-pi/slice_count ...
                                powder_grid.betas(crystal_idx) ...
                                powder_grid.gammas(crystal_idx)];
        liouv_stack=rotor_stack(spin_system,parameters,'labframe');

        % Propagate towards decreasing rotor phase
        rho=parameters.rho0;
        for slice_idx=1:slice_count
            rotor_idx=mod(1-slice_idx,slice_count)+1;
            rho=expm(-1i*full(liouv_stack{rotor_idx}+parameters.rf_amp*...
                                parameters.rf_op)*...
                     (parameters.duration/slice_count))*rho;
        end

        % Accumulate the rotor-phase signal
        crystal_signal=crystal_signal+parameters.coil'*rho/phase_count;
    end

    % Apply the crystallite weight
    signal=signal+powder_grid.weights(crystal_idx)*...
                  crystal_signal/parameters.ref_norm;
end

end

% Add a constant RF generator across FP rotor phase collocation points
function signal=fp_powder_signal(~,parameters,generator,~,~)

% Add the RF operator at every rotor collocation point
generator=generator+parameters.rf_amp*kron(speye(parameters.spc_dim),...
                           parameters.rf_op);

% Propagate the uniform rotor-phase state for one period
rho=expm(-1i*full(generator)*parameters.duration)*parameters.rho0;

% Detect and normalise the signal
signal=(parameters.coil'*rho)/parameters.ref_norm;

end


