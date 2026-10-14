% Snapshot POD for repeated spectra of a coupled two-proton system.
% Five explicitly sampled offset trajectories train one Euclidean basis;
% interleaved held-out offsets compare full Spinach and projected dynamics.
% A rank-three truncation and the full four-dimensional coherence subspace
% are both reported: snapshot energy is not an observable error bound.
% Training and projection costs must be amortised separately from evaluation.
%
% ilya.kuprov@weizmann.ac.il

function snapshot_pod_nmr()

% Specify a strongly coupled proton pair with transverse relaxation
sys.magnet=14.1; sys.isotopes={'1H','1H'}; sys.output='hush';
sys.parallel={'local',1}; sys.disable={'hygiene'};
inter.zeeman.scalar={0 0.25}; inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=35;
inter.relaxation={'t1_t2'}; inter.r1_rates={0.7 0.8};
inter.r2_rates={2 2.5}; inter.equilibrium='zero'; inter.rlx_keep='labframe';
bas.formalism='sphten-liouv'; bas.approximation={'none'};

% Build the complete Liouville system and an explicit affine offset term
spin_system=create(sys,inter); spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');
spin_system.tols.prop_chop=eps; spin_system.tols.liouv_zero=eps;
spin_system.tols.zte_tol=eps;
H=hamiltonian(spin_system); R=relaxation(spin_system);
terms={H+1i*R;2*pi*operator(spin_system,'Lz',1)};
rho=state(spin_system,'L+',1); coil=rho;

% Specify the training coverage and independent acquisition samples
train_offsets=[-30 -15 0 15 30]; test_offsets=[-22.5 -7.5 7.5 22.5];
timestep=1/1000; nsteps=255;
options=struct('method','svd','oversampling',0,'power_iters',0,'seed',0);
blocks=cell(numel(train_offsets),1);
for n=1:numel(train_offsets)
    L=terms{1}+train_offsets(n)*terms{2};
    blocks{n}=full(evolution(spin_system,L,[],rho,timestep,nsteps,'trajectory'));
end

% Train explicit ranks and measure held-out FIDs and identically processed spectra
window=cos(linspace(0,pi/2,nsteps+1)).^2;
for target_rank=[3 4]
    prepared=pod_basis(blocks,ones(size(blocks)),target_rank,options);
    model=pod_project(prepared.basis,terms,full(rho),full(coil));
    for n=1:numel(test_offsets)
        L=terms{1}+test_offsets(n)*terms{2};
        full_fid=evolution(spin_system,L,coil,rho,timestep,nsteps,'observable').';
        reduced_fid=pod_evolve(model,[1;test_offsets(n)],timestep,nsteps);
        full_spec=fftshift(fft(full_fid.*window,1024));
        reduced_spec=fftshift(fft(reduced_fid.*window,1024));
        peak=max(abs(full_spec)); delta=reduced_spec-full_spec;
        fprintf('rank=%d offset=%+.1f Hz fid_max=%.3e peak_max=%.3e peak_rms=%.3e\n',...
                target_rank,test_offsets(n),max(abs(reduced_fid-full_fid)),...
                max(abs(delta))/peak,sqrt(mean(abs(delta).^2))/peak);
    end
    fprintf('rank=%d supplied_snapshot_residual=%.3e\n',target_rank,prepared.residual_rel);
end

% Display the last independently processed full and reduced spectra on one scale
frequency=(-512:511)/(1024*timestep);
kfigure(); plot(frequency,real(full_spec),'k-',frequency,real(reduced_spec),'r--');
kxlabel('frequency, Hz'); kylabel('real spectral amplitude'); kgrid;
klegend({'Full Spinach','Snapshot POD'},'Location','Best');

end

