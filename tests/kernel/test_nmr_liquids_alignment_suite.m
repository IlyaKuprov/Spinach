% Tests compact literature-alignment probes for liquid-state NMR
% pulse sequences. Syntax:
%
%              result=test_nmr_liquids_alignment_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test covers compact gCOSY, HMBC, HSQC, and TOCSY paths that
% were updated during the nmr_liquids literature-alignment pass.
%
% ilya.kuprov@weizmann.ac.il

function result=test_nmr_liquids_alignment_suite()

% Announce the test target
fprintf('TESTING: NMR liquids literature-alignment probes\n');

% State the test target
result=new_test_result('kernel/nmr_liquids_alignment_suite',...
                       'NMR liquids literature-alignment probes',...
                       'Updated liquid-state pulse sequence paths must produce finite non-zero compact FIDs.');

% Build a compact homonuclear system for gCOSY
sys.magnet=9.4;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=8;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Set compact gCOSY parameters
parameters.sweep=100;
parameters.npoints=[4 4];
parameters.spins={'1H'};
parameters.angle=pi/2;
parameters.g_amp=3;
parameters.g_dur=1e-3;
parameters.g_stab_del=0;
parameters.s_len=1.5;
parameters.pathway='P';

% Run the P-type gCOSY pathway
fid_p=liquid(spin_system,@gcosy,parameters,'nmr');

% Run the N-type gCOSY pathway
parameters.pathway='N';
fid_n=liquid(spin_system,@gcosy,parameters,'nmr');

% Run the phase-sensitive gCOSY pathway pair
parameters.pathway='P+N';
fid_pn=liquid(spin_system,@gcosy,parameters,'nmr');

% Check gCOSY outputs
result=test_close(result,'gCOSY P finite',all(isfinite(fid_p(:))),true,0,0,...
                  'the P-type gradient pathway must produce finite data');
result=test_close(result,'gCOSY N finite',all(isfinite(fid_n(:))),true,0,0,...
                  'the N-type gradient pathway must produce finite data');
result=test_close(result,'gCOSY pathways non-zero',...
                  (norm(fid_p(:))>0)&&(norm(fid_n(:))>0),true,0,0,...
                  'both gradient-selected pathways must retain observable signal');
result=test_close(result,'gCOSY pathways distinct',...
                  norm(fid_p(:)-fid_n(:))>0,true,0,0,...
                  'P-type and N-type gradient pathways must not be identical');
result=test_close(result,'gCOSY P+N structure',...
                  isstruct(fid_pn)&&...
                  isfield(fid_pn,'pos')&&isfield(fid_pn,'neg'),true,0,0,...
                  'P+N mode must return both gradient-selected pathway components');
result=test_close(result,'gCOSY P+N P branch',fid_pn.pos,fid_p,1e-10,1e-10,...
                  'P+N mode P-type branch must match direct P acquisition');
result=test_close(result,'gCOSY P+N N branch',fid_pn.neg,fid_n,1e-10,1e-10,...
                  'P+N mode N-type branch must match direct N acquisition');

% Recombine the gCOSY pathways for phase-sensitive processing
f1_pos=fftshift(fft(fid_pn.pos,parameters.npoints(2),1),1);
f1_neg=fftshift(fft(fid_pn.neg,parameters.npoints(2),1),1);
fid_states=f1_pos+conj(f1_neg);
spec_pn=fftshift(fft(fid_states,parameters.npoints(1),2),2);

% Check the recombined gCOSY spectrum
result=test_close(result,'gCOSY P+N recombined finite',...
                  all(isfinite(spec_pn(:)))&&(norm(real(spec_pn(:)))>0),true,0,0,...
                  'P+N pathway recombination must produce finite non-zero data');

% Build a compact heteronuclear system for HMBC and HSQC
sys.isotopes={'1H','15N'};
inter.zeeman.scalar={0 0};
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=90;
spin_system=test_spin_system(sys,inter,bas);

% Set compact HMBC parameters
parameters=struct();
parameters.sweep=[100 100];
parameters.npoints=[4 4];
parameters.spins={'15N','1H'};
parameters.J=90;
parameters.delta_b=5e-3;

% Run HMBC on a non-carbon heteronucleus
fid_hmbc=liquid(spin_system,@hmbc,parameters,'nmr');

% Check HMBC output
result=test_close(result,'HMBC non-carbon finite',...
                  all(isfinite(fid_hmbc(:)))&&(norm(fid_hmbc(:))>0),true,0,0,...
                  'HMBC coherence selection must follow parameters.spins{1}');

% Set compact HSQC parameters
parameters.decouple_f1={'1H'};
parameters.decouple_f2={'15N'};

% Run HSQC
fid_hsqc=liquid(spin_system,@hsqc,parameters,'nmr');

% Check HSQC outputs
result=test_close(result,'HSQC positive pathway finite',...
                  all(isfinite(fid_hsqc.pos(:)))&&(norm(fid_hsqc.pos(:))>0),true,0,0,...
                  'HSQC positive pathway must produce finite observable data');
result=test_close(result,'HSQC negative pathway finite',...
                  all(isfinite(fid_hsqc.neg(:)))&&(norm(fid_hsqc.neg(:))>0),true,0,0,...
                  'HSQC negative pathway must produce finite observable data');

% Run CT-HSQC
fid_ct_hsqc=liquid(spin_system,@ct_hsqc,parameters,'nmr');

% Check CT-HSQC outputs
result=test_close(result,'CT-HSQC positive pathway finite',...
                  all(isfinite(fid_ct_hsqc.pos(:)))&&...
                  isequal(size(fid_ct_hsqc.pos),[parameters.npoints(2) parameters.npoints(1)])&&...
                  (norm(fid_ct_hsqc.pos(:))>0),true,0,0,...
                  'CT-HSQC positive pathway must produce finite observable data');
result=test_close(result,'CT-HSQC negative pathway finite',...
                  all(isfinite(fid_ct_hsqc.neg(:)))&&...
                  isequal(size(fid_ct_hsqc.neg),[parameters.npoints(2) parameters.npoints(1)])&&...
                  (norm(fid_ct_hsqc.neg(:))>0),true,0,0,...
                  'CT-HSQC negative pathway must produce finite observable data');

% Set compact NOESY-HSQC parameters
parameters=struct();
parameters.sweep=[100 100 100];
parameters.npoints=[2 2 2];
parameters.spins={'1H','15N','1H'};
parameters.J=90;
parameters.tmix=1e-3;

% Run NOESY-HSQC through the decoupled mixing path
fid_noesyhsqc=liquid(spin_system,@noesyhsqc,parameters,'nmr');

% Check NOESY-HSQC outputs
result=test_close(result,'NOESY-HSQC finite',...
                  all(isfinite(fid_noesyhsqc.pos_pos(:)))&&...
                  all(isfinite(fid_noesyhsqc.pos_neg(:)))&&...
                  all(isfinite(fid_noesyhsqc.neg_pos(:)))&&...
                  all(isfinite(fid_noesyhsqc.neg_neg(:))),true,0,0,...
                  'NOESY-HSQC decoupled mixing path must produce finite data');

% Build a compact direct TOCSY probe
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.coupling.scalar=cell(1,1);
spin_system=test_spin_system(sys,inter,bas);
rho0=state(spin_system,'Lz','1H','cheap');
dim=size(rho0,1);
H=sparse(dim,dim);
R=-speye(dim);
K=sparse(dim,dim);

% Set compact TOCSY parameters
parameters=struct();
parameters.sweep=[100 100];
parameters.npoints=[4 4];
parameters.spins={'1H'};
parameters.tmix=0;
parameters.lamp=1000;
parameters.rho0=rho0;

% Run TOCSY without relaxation time
fid_ref=tocsy(spin_system,parameters,H,R,K);

% Run TOCSY with relaxation time
parameters.tmix=0.2;
fid_rel=tocsy(spin_system,parameters,H,R,K);

% Check that relaxation acts during the spin-lock interval
result=test_close(result,'TOCSY spin-lock relaxation attenuation',...
                  norm([fid_rel.cos(:); fid_rel.sin(:)])<...
                  norm([fid_ref.cos(:); fid_ref.sin(:)]),true,0,0,...
                  'TOCSY mixing must include relaxation during the spin-lock interval');

end
