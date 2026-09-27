% Tests cache identity across Hilbert-to-Liouville conversion. Syntax:
%
%                      result=test_sim2liouv_cache()
%
% Outputs:
%
%     result - regression checks for both cache insertion orders,
%              complex operators, soft-pulse acquisition, and no-op paths
%
% talos@spindynamics.org

function result=test_sim2liouv_cache()

% State the cache-conversion target
result=new_test_result('kernel/sim2liouv_cache','Converted basis cache identity',...
                       'converted operators and Hamiltonians must use the converted basis cache identity.');

% Build small systems with a real one-worker cache store
sys.magnet=1; sys.isotopes={'1H'}; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.parprops={};
inter.zeeman.scalar={1}; bas.formalism='zeeman-hilb'; bas.approximation='none';
flags={{'op_cache'},{'ham_cache'},{'op_cache','ham_cache'}};
nuclei={'1H','13C'};
for flag_idx=1:numel(flags)
    sys.enable=flags{flag_idx};
    for ordering=1:2

        % Separate the fixtures without altering physical cache keys
        sys.isotopes=nuclei(ordering); inter.zeeman.scalar={flag_idx};
        spin_h=assume(basis(create(sys,inter),bas),'nmr');
        ref_h=spin_h; ref_h.sys.enable={}; H=hamiltonian(ref_h);
        [spin_l,~,K]=sim2liouv(spin_h,struct(),H,[],[]);
        ref_l=spin_l; ref_l.sys.enable={};
        result=test_close(result,'canonical converted hash',...
                          strcmp(spin_l.bas.basis_hash,md5_hash(spin_l.bas.basis)),true,0,0,...
                          'the converted basis must have its canonical identity');

        % Warm both representations in the requested order
        if ordering==1
            operator(spin_h,'L+',1); hamiltonian(spin_h);
            op_obs=operator(spin_l,'L+',1); h_obs=hamiltonian(spin_l);
            op_ref=operator(ref_l,'L+',1); h_ref=K;
        else
            operator(spin_l,'L+',1); hamiltonian(spin_l);
            op_obs=operator(spin_h,'L+',1); h_obs=hamiltonian(spin_h);
            op_ref=operator(ref_h,'L+',1); h_ref=H;
        end
        result=test_close(result,'cached raising operator',op_obs,op_ref,1e-12,1e-12,...
                          'a non-Hermitian operator must not collide across representations');
        result=test_close(result,'cached drift',h_obs,h_ref,1e-12,1e-12,...
                          'Hamiltonians must retain their requested representation');
    end
end

% Keep existing cache metadata correct while caching is temporarily disabled
spin_h.sys.enable={};
[spin_l,~,~,~,~]=sim2liouv(spin_h,struct(),[],[],[]);
result=test_close(result,'disabled cache hash',...
                  strcmp(spin_l.bas.basis_hash,md5_hash(spin_l.bas.basis)),true,0,0,...
                  'existing metadata must remain valid if caching is enabled again');

% Preserve objects that never requested a cache
sys.enable={}; spin_h=basis(create(sys,inter),bas);
[spin_l,~,~,~,~]=sim2liouv(spin_h,struct(),[],[],[]);
result=test_close(result,'absent cache metadata',isfield(spin_l.bas,'basis_hash'),false,0,0,...
                  'conversion without caching must not introduce cache metadata');

% Return every input unchanged for the three no-op formalisms
forms={'zeeman-liouv','sphten-liouv','zeeman-wavef'};
for form_idx=1:numel(forms)
    bas.formalism=forms{form_idx}; sys.enable={'op_cache','ham_cache'};
    spin_system=basis(create(sys,inter),bas); parameters.tag='unchanged';
    H=operator(spin_system,'Lx',1)+0.37*operator(spin_system,'Ly',1);
    R=sparse(size(H,1),size(H,2)); K=R;
    [observed,params_out,h_out,r_out,k_out]=sim2liouv(spin_system,parameters,H,R,K);
    result=test_close(result,'no-op inputs',...
                      isequal(observed,spin_system)&&isequal(params_out,parameters)&&...
                      isequal(h_out,H)&&isequal(r_out,R)&&isequal(k_out,K),true,0,0,...
                      'non-Hilbert formalisms must be returned unchanged');
end

% Acquire a complex FID after a noncommuting phase-shifted soft pulse
sys.magnet=1; sys.isotopes={'1H'}; sys.enable={'op_cache','ham_cache'};
inter.zeeman.scalar={0}; bas.formalism='zeeman-hilb';
spin_system=basis(create(sys,inter),bas);
parameters=struct(); parameters.spins={'1H'};
parameters.rho0=state(spin_system,'Lz','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.pulse_frq=17; parameters.pulse_phi=0.37;
parameters.pulse_pwr=2*pi*100; parameters.pulse_dur=0.0025;
parameters.pulse_rnk=2; parameters.offset=31;
parameters.sweep=1000; parameters.npoints=8; parameters.method='expm';
fid_cached=liquid(spin_system,@sp_acquire,parameters,'nmr');
spin_system.sys.enable={};
fid_ref=liquid(spin_system,@sp_acquire,parameters,'nmr');
result=test_close(result,'cached soft-pulse FID',fid_cached,fid_ref,1e-11,1e-11,...
                  'cache use must not change the acquired signal');
result=test_close(result,'nonzero complex signal',...
                  norm(fid_ref)>0.1&&norm(imag(fid_ref))>0.1,true,0,0,...
                  'the acquisition must exercise populated complex dynamics');

end


