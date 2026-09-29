% Compares powder-averaged sliced MAS with Fokker-Planck MAS dynamics.
% Syntax:
%
%                        mas_fplanck_powder()
%
% No inputs. Checks rotor rank, slices, phase, and powder grid convergence.
% The same weighted Lebedev crystallite grid is used by both routes.
%
% Calculation time: minutes
%
% talos@spindynamics.org

function mas_fplanck_powder()

% Set a spin with anisotropic shielding and continuous transverse RF
sys.magnet=9.4;
sys.isotopes={'13C'};
sys.parallel={'processes',1};
sys.parprops={};
inter.zeeman.eigs={[-120 -25 145]};
inter.zeeman.euler={[0.4 0.7 0.2]};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system.sys.output='hush';
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
% The complete rotor stack spans exactly one rotor period
parameters.duration=1/parameters.rate;
parameters.serial=true;
parameters.verbose=0;

% Refine the Fokker-Planck rotor rank on the fixed powder grid
ranks=[6 8];
fp_sig=zeros(size(ranks));
for n=1:numel(ranks)
    parameters.max_rank=ranks(n);
    fp_sig(n)=singlerot(spin_system,@fp_powder_signal,...
                       parameters,'labframe');
    fprintf('Powder FP rank %d: %.9g%+.9gi\n',ranks(n),...
            real(fp_sig(n)),imag(fp_sig(n)));
end

% Refine phase quadrature and midpoint slices independently
counts=[33 65]; phase_counts=[7 13];
sl_sig=zeros(numel(counts),numel(phase_counts));
for n=1:numel(counts)
    for k=1:numel(phase_counts)
        sl_sig(n,k)=sliced_powder_signal(spin_system,parameters,...
                                        counts(n),phase_counts(k));
        fprintf('Powder slices %d phases %d: %.9g%+.9gi\n',...
                counts(n),phase_counts(k),...
                real(sl_sig(n,k)),imag(sl_sig(n,k)));
    end
end

% Check the powder quadrature against a finer Lebedev grid
parameters.grid='leb_2ang_rank_11';
fp_fine=singlerot(spin_system,@fp_powder_signal,...
                  parameters,'labframe');
grid_step=abs(fp_fine-fp_sig(end));

% Compare normalised signals and each independent refinement
fp_step=abs(fp_sig(end)-fp_sig(end-1));
sl_step=abs(sl_sig(end,end)-sl_sig(end-1,end));
phase_step=abs(sl_sig(end,end)-sl_sig(end,end-1));
route_gap=abs(fp_sig(end)-sl_sig(end,end));
fprintf('Powder CSA: FP step %.6g, slice step %.6g, ',fp_step,sl_step);
fprintf('phase step %.6g, grid step %.6g, gap %.6g\n',...
        phase_step,grid_step,route_gap);
target=0.002;
assert(all(isfinite([fp_step sl_step phase_step grid_step ...
                     route_gap]))&&...
       fp_step<target/2&&sl_step<target/2&&...
       phase_step<target/2&&grid_step<target&&route_gap<target,...
       'Powder-averaged MAS routes did not converge.');
fprintf('MAS_FPLANCK_POWDER_SUCCESS\n');

end

% Build the sliced propagator for each crystallite and rotor start phase
function signal=sliced_powder_signal(spin_system,parameters,count,nphases)

grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',...
                   parameters.grid),'alphas','betas','gammas','weights');
parameters.grid='single_crystal';
parameters.masframe='rotor';
parameters.max_rank=(count-1)/2;
signal=0;
for q=1:numel(grid.weights)
    orient_sig=0;
    for phase=1:nphases
        rotor_phase=2*pi*(phase-1)/nphases;
        parameters.orientation=[grid.alphas(q)+rotor_phase-pi/count ...
                                grid.betas(q) grid.gammas(q)];
        L=rotor_stack(spin_system,parameters,'labframe');
        rho=parameters.rho0;
        for n=1:count
            idx=mod(1-n,count)+1;
            rho=expm(-1i*full(L{idx}+parameters.rf_amp*...
                                parameters.rf_op)*...
                     (parameters.duration/count))*rho;
        end
        orient_sig=orient_sig+parameters.coil'*rho/nphases;
    end
    signal=signal+grid.weights(q)*orient_sig/parameters.ref_norm;
end

end

% Add a constant RF generator across FP rotor phase collocation points
function signal=fp_powder_signal(~,parameters,G,~,~)

G=G+parameters.rf_amp*kron(speye(parameters.spc_dim),...
                           parameters.rf_op);
rho=expm(-1i*full(G)*parameters.duration)*parameters.rho0;
signal=(parameters.coil'*rho)/parameters.ref_norm;

end


