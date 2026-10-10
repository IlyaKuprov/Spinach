% Tests unweighted sequence probes and weighted PRESS preparation. Syntax:
%
%                      result=test_cwdm_sequences()
%
% Outputs:
%
%    result - concentration scaling and zero-population checks
%
% Normalised ENDOR and relaxation probes must not vanish at zero population.
% PRESS phantom amplitudes are linear, not quadratic, in concentration.
% Small homogeneous grids exercise all four migrated imaging sequences.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_sequences()

% Announce the sequence contract
fprintf('TESTING: CWDM sequence probes\n');
result=new_test_result('kernel/cwdm_sequences','CWDM sequence probes',...
                      'Receivers and normalised probes do not carry concentrations.');

% Build a proton with known longitudinal and transverse rates
sys.magnet=1; sys.isotopes={'1H'};
inter.zeeman.scalar={0}; inter.relaxation={'t1_t2'};
inter.r1_rates={2}; inter.r2_rates={5}; inter.equilibrium='zero';
inter.rlx_keep='diagonal';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
s=test_spin_system(sys,inter,bas); empty=s; empty.chem.concs=0;
[r1,r2]=relaxan(empty);
result=test_close(result,'empty longitudinal rate',r1,2,1e-12,0,...
                  'the longitudinal Rayleigh quotient remains defined at zero population');
result=test_close(result,'empty transverse rate',r2,5,1e-12,0,...
                  'the transverse Rayleigh quotient remains defined at zero population');

% Compare PRESS profiles at unit, fractional, and zero concentrations
for ndim=2:3
    parameters.spins={'1H'}; parameters.npts=2*ones(1,ndim);
    parameters.ss_grad_amp=zeros(1,ndim);
    parameters.rf_frq_list=repmat({0},1,ndim);
    parameters.rf_amp_list=repmat({pi/2},1,ndim);
    parameters.rf_dur_list=repmat({1},1,ndim);
    parameters.rf_phi=repmat({0},1,ndim);
    parameters.max_rank=repmat({2},1,ndim);
    Z=sparse(4*prod(parameters.npts),4*prod(parameters.npts));
    H=polyadic({{opium(prod(parameters.npts),1),operator(s,'Lz','1H')}});
    G=repmat({Z},1,ndim); sequence=str2func(sprintf('press_voxel_%dd',ndim));
    reference=sequence(s,parameters,H,Z,Z,G,Z);
    for conc=[0.3 0]
        local=s; local.chem.concs=conc;
        actual=sequence(local,parameters,H,Z,Z,G,Z);
        result=test_close(result,sprintf('PRESS %dD concentration %g',ndim,conc),...
                          actual,conc*reference,1e-12,1e-12,...
                          'initial density is weighted once and the final probe is unweighted');
    end
    result=test_true(result,sprintf('PRESS %dD nonzero reference',ndim),...
                     norm(reference(:))>0,'the concentration comparison exercises a nonzero signal');
end

% Exercise both three-dimensional slice displays and acquisitions
parameters=struct(); parameters.npts=[2 2 2]; parameters.dims=[1 1 1];
parameters.image_size=[3 3]; parameters.ss_grad_amp=0;
parameters.pe_grad_amp=0; parameters.ro_grad_amp=0;
parameters.pe_grad_dur=0.01; parameters.ro_grad_dur=0.01;
parameters.t_echo=0.01; parameters.rf_frq_list=0;
parameters.rf_amp_list=pi/2; parameters.rf_dur_list=1; parameters.rf_phi=0;
parameters.rho0=kron(ones(8,1),state(s,'Lz','1H'));
parameters.coil=kron(ones(8,1),coil_state(s,'L+','1H','exact'));
Z=sparse(32,32); G={Z,Z,Z};
for sequence={@epi_3d,@phase_enc_3d}
    reference=sequence{1}(s,parameters,Z,Z,Z,G,Z);
    actual=sequence{1}(empty,parameters,Z,Z,Z,G,Z);
    result=test_close(result,[func2str(sequence{1}) ' empty geometry'],actual,reference,...
                      1e-12,1e-12,'caller-supplied density and receiver do not change with stored concentration');
end
close all;

% Build a hyperfine-coupled electron and proton for normalised CW ENDOR
sys.isotopes={'E','1H'}; inter=struct();
inter.zeeman.scalar={2.0023,0}; inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=1e6;
s=test_spin_system(sys,inter,bas); empty=s; empty.chem.concs=0;
parameters=struct('sweep',1e6,'npoints',4); Z=sparse(16,16);
reference=endor_cw(s,parameters,Z,Z,Z);
actual=endor_cw(empty,parameters,Z,Z,Z);
result=test_close(result,'normalised ENDOR zero population',actual,reference,1e-12,1e-12,...
                  'normalisation uses an unweighted nuclear operator sum');
result=test_true(result,'ENDOR nonzero reference',norm(reference)>0,...
                 'the normalised spectrum is nonzero');

% Hold level populations fixed while comparing microwave transition geometry
sys.isotopes={'E'}; inter=struct('zeeman',struct('matrix',{{2.0023*eye(3)}}));
bas.formalism='zeeman-hilb'; s=test_spin_system(sys,inter,bas);
parameters=struct('spins',{{'E'}},'grid','icos_2ang_12pts',...
                  'mw_freq',9.5e9,'fwhm',2e-3,'window',[0.33 0.35],...
                  'npoints',9,'tm_tol',0,'rspt_order',Inf,'int_tol',1e9);
parameters.rho0=-coil_state(s,'Lz','E','exact');
reference=fieldsweep(s,parameters); s.chem.concs=0;
actual=fieldsweep(s,parameters);
result=test_close(result,'microwave operator zero population',actual,reference,1e-12,1e-12,...
                  'the transition operator is independent of stored concentration');
result=test_true(result,'field sweep nonzero reference',norm(reference)>0,...
                 'the isotropic electron produces a nonzero field-swept signal');
fprintf('CWDM_SEQUENCE_PATHS_COMPLETE functions=7 substitutions=8\n');

end


