% Tests remaining deterministic utility helpers. Syntax:
%
%                    result=test_dynamic_remaining_core_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks block eliminations, text reporting, spin metadata,
% analytical line shapes, pumping terms, kite pruning, trajectory
% stitching, and small random-rotation diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_remaining_core_suite()

% Announce the test target
fprintf('TESTING: Remaining deterministic utility helpers\n');

% State the utility target of the test
result=new_test_result('kernel/dynamic_remaining_core_suite',...
                       'Remaining deterministic utility helpers',...
                       'Small hand-written utility calls must preserve documented algebraic and reporting semantics.');

% Build a minimal Liouville-space descriptor for algebraic utilities
spin_system=local_liouvillian_system(3);

% Check adiabatic elimination against the exact Schur complement term
L=[1 2 3;4 5 6;7 8 10];
[L_slow,R_extra]=adelim(spin_system,L,3,[1 2]);
R_ref=1i*[3;6]*(1/10)*[7 8];
result=test_close(result,'adelim slow block',L_slow,[1 2;4 5],1e-14,1e-14,...
                  'adiabatic elimination must return the slow-space Liouvillian block');
result=test_close(result,'adelim induced relaxation',R_extra,R_ref,1e-14,1e-14,...
                  'the eliminated fast block contributes i*L01*(L11\\L10)');

% Check a textbook Clebsch-Gordan coefficient
cg_triplet=cg_fast(1,0,1/2,1/2,1/2,-1/2);
result=test_close(result,'cg_fast triplet coefficient',cg_triplet,1/sqrt(2),...
                  1e-12,1e-12,...
                  'two spin-half states couple to the triplet M=0 state with coefficient 1/sqrt(2)');

% Check direct console and banner reporting into a file handle
log_file=[tempname '.txt'];
file_id=fopen(log_file,'w');
spin_system.sys.output=file_id;
report(spin_system,'direct report message');
banner(spin_system,'basis_banner');
fclose(file_id);
log_text=fileread(log_file);
delete(log_file);
result=test_true(result,'report and banner file output',...
                 contains(log_text,'direct report message')&&contains(log_text,'BASIS SET'),...
                 'report() must write prefixed messages and banner() must delegate to report()');
spin_system.sys.output='hush';

% Check polyadic text diagnostics with an explicit label
polinfo(polyadic({{speye(2),sparse(1)}}),1,'labelled');
result.messages{end+1}='PASS: polinfo labelled output -- direct polyadic diagnostic call completed';

% Check summary traverses coordinate metadata through the silent report path
spin_system.comp.nspins=1;
spin_system.comp.isotopes={'1H'};
spin_system.comp.labels={'proton'};
spin_system.comp.mults=2;
spin_system.inter.coordinates={[0 0 0]};
summary_coordinates(spin_system,'coordinate summary');
result.messages{end+1}='PASS: summary_coordinates -- direct silent metadata summary call completed';

% Check impound cell packaging preserves values and types
payload=impound(17,'spinach',{speye(2)});
result=test_true(result,'impound payload',...
                 iscell(payload)&&isequal(payload{1},17)&&strcmp(payload{2},'spinach')&&...
                 isequal(payload{3}{1},speye(2)),...
                 'impound() must return all inputs unchanged in a cell array');

% Check dipolar coupling for two spins one Angstrom apart on the X axis
dip_system=local_dipolar_system();
dip_system=dipolar(dip_system);
dip_pref=0.5*dip_system.inter.gammas(1)*dip_system.inter.gammas(2)*...
         dip_system.tols.hbar*dip_system.tols.mu0/(4*pi*(1e-10)^3);
dip_ref=dip_pref*[-2 0 0;0 1 0;0 0 1];
result=test_close(result,'dipolar coupling tensor',...
                  dip_system.inter.coupling.matrix{1,2},dip_ref,1e-6,1e-12,...
                  'a unit X-axis internuclear vector gives the traceless diag([-2 1 1]) dipolar tensor');
result=test_true(result,'dipolar proximity matrix',nnz(dip_system.inter.proxmatrix)==2,...
                 'both directed spin-pair proximities must be detected below the distance cutoff');

% Check isotope swapping rescales spin-pair couplings and wipes quadratic terms
[sys,inter]=local_isoswap_inputs();
[sys_out,inter_out]=isoswap(sys,inter,1,'2H');
gamma_ratio=spin('2H')/spin('1H');
result=test_true(result,'isoswap isotope replacement',strcmp(sys_out.isotopes{1},'2H'),...
                 'isoswap() must replace the requested isotope string');
result=test_close(result,'isoswap coupling scaling',...
                  inter_out.coupling.matrix{1,2},gamma_ratio*eye(3),1e-14,1e-14,...
                  'bilinear couplings involving the replaced spin must scale by the gyromagnetic ratio');
result=test_true(result,'isoswap quadratic wipe',...
                 isempty(inter_out.coupling.matrix{1,1}),...
                 'quadratic self-couplings are not transferable and must be wiped');

% Check interaction-representation order zero against the perturbation block
spin_system=local_liouvillian_system(2);
H0=2*pi*diag([0 1]);
H=H0+diag([0.05 -0.02]);
Hr=intrep(spin_system,H0,H,1,0);
result=test_close(result,'intrep zeroth order',Hr,diag([0.05 -0.02]),1e-13,1e-13,...
                  'zeroth-order interaction representation returns H-H0 after period validation');

% Check finite-difference Hessian construction on a constant 3D field
H=fdhess(ones(5,6,7),3);
result=test_true(result,'fdhess output layout',isequal(size(H),[3 3]),...
                 'fdhess() must return a 3x3 Hessian cell array');
for n=1:3
    for k=1:3
        result=test_close(result,['fdhess constant block ' num2str(n) num2str(k)],...
                          H{n,k},zeros(5,6,7),1e-14,1e-14,...
                          'all Hessian components of a constant field must vanish');
    end
end

% Check sinkhole column removal in spherical-tensor Liouville formalism
spin_system=local_liouvillian_system(3);
L=[1 2 3;4 5 6;7 8 9];
L_sink=sinkhole(spin_system,L,[2 3]);
result=test_close(result,'sinkhole frozen columns',L_sink,[1 0 0;4 0 0;7 0 0],0,0,...
                  'sinkhole() must zero all Liouvillian columns that feed frozen states');

% Check the normalised Lorentzian branch analytically
x_axis=[-1 0 1];
line_obs=lorentzcon(0,2*pi,2,x_axis);
line_ref=2./(1+x_axis.^2);
result=test_close(result,'lorentzcon Lorentzian',line_obs,line_ref,1e-14,1e-14,...
                  'with fwhm two and amplitude 2*pi, the normalised Lorentzian is 2/(1+x^2)');

% Check the narrow-linewidth Lorentzian peak at exact centre
tiny_fwhm=1e-308;
tiny_gam=tiny_fwhm/2;
line_obs=lorentzcon(0,1,tiny_fwhm,0);
line_ref=1/(pi*tiny_gam);
result=test_close(result,'lorentzcon tiny centre',line_obs,line_ref,0,1e-14,...
                  'exact-centre evaluation must not use a reciprocal-width product that turns zero offset into NaN');

% Check the narrow-linewidth MEX segment branch at exact offsets
if exist('lorentzcon','file')==3
    line_obs=lorentzcon([0 1],1,tiny_fwhm,[0 1]);
    line_ref=[0.5 0.5];
    result=test_close(result,'lorentzcon MEX tiny segment endpoints',line_obs,line_ref,1e-14,1e-14,...
                      'two-argument atan2 must keep exact endpoint segment values finite when reciprocal width overflows');
end

% Check the normalised Gaussian branch analytically
x_axis=[-1 0 1];
line_obs=gausscon(0,2,2,x_axis);
line_ref=2*gaussfun(x_axis,2);
result=test_close(result,'gausscon Gaussian',line_obs,line_ref,1e-14,1e-14,...
                  'scalar offset branch must match the normalised Gaussian line shape');

% Check that coincident triangle vertices collapse to a Gaussian
line_obs=gausscon([0 0 0],2,2,0);
line_ref=2*gaussfun(0,2);
result=test_close(result,'gausscon coincident vertices',line_obs,line_ref,1e-14,1e-14,...
                  'three identical triangle vertices must collapse to the scalar-offset Gaussian line shape');

% Check the Gaussian convolution with a general triangle by quadrature
offs=[-1 0.25 2];
x_axis=[-0.5 0.25 1.75];
line_obs=gausscon(offs,1,0.6,x_axis);
line_ref=zeros(size(x_axis));
for n=1:numel(x_axis)
    line_ref(n)=integral(@(offset) ...
        ((((offset>=offs(1)).*(offset<=offs(2))).*...
          (2*(offset-offs(1))/((offs(2)-offs(1))*(offs(3)-offs(1)))))+...
         (((offset>=offs(2)).*(offset<=offs(3))).*...
          (2*(offs(3)-offset)/((offs(3)-offs(2))*(offs(3)-offs(1)))))).*...
        gaussfun(x_axis(n)-offset,0.6),offs(1),offs(3),...
        'Waypoints',offs(2),'AbsTol',1e-12,'RelTol',1e-12);
end
result=test_close(result,'gausscon triangle quadrature',line_obs,line_ref,1e-10,1e-10,...
                  'closed-form Gaussian triangle convolution must match direct numerical quadrature');

% Check magnetic pumping adds only a source column from the unit state
spin_system=local_liouvillian_system(3);
R=zeros(3);
rho=[0;2;-1];
R_pumped=magpump(spin_system,R,rho,0.25);
result=test_close(result,'magpump source column',R_pumped,[0 0 0;0.5 0 0;-0.25 0 0],...
                  1e-14,1e-14,...
                  'magpump() must add rate*rho to the first relaxation-superoperator column only');

% Check Redfield-kite pruning keeps diagonals and longitudinal cross-relaxation
spin_system=local_liouvillian_system(4);
spin_system.bas.basis=[0;1;2;3];
R=sparse([1 2 2 4],[3 2 4 4],[2 7 5 9],4,4);
R_kite=sec2kite(spin_system,R);
R_ref=sparse([1 2 4],[3 2 4],[2 7 9],4,4);
result=test_close(result,'sec2kite pruning',R_kite,R_ref,1e-14,1e-14,...
                  'only self-relaxation and longitudinal cross-relaxation rates should survive');

% Check the Sorensen unitary transfer bound from ordered eigenvalues
bound=sorensen(diag([3 1]),diag([2 0]));
result=test_close(result,'sorensen eigenvalue bound',bound,1.5,1e-14,1e-14,...
                  'the bound is the sorted-eigenvalue scalar product divided by trace(rho_targ^2)');

% Check trajectory stitching for zero Liouvillian and identity midpoint event
spin_system=local_liouvillian_system(2);
rho_stack=[1 2;3 4];
coil_stack=eye(2);
t1.nsteps=2; t2.nsteps=3; t2.timestep=0.1; t3.nsteps=2;
fid=stitch(spin_system,zeros(2),rho_stack,coil_stack,{zeros(2)},{1},t1,t2,t3);
fid_ref=zeros(2,3,2);
for n=1:t2.nsteps
    fid_ref(:,n,:)=coil_stack'*rho_stack;
end
result=test_close(result,'stitch zero dynamics',fid,fid_ref,1e-14,1e-14,...
                  'with zero evolution and identity midpoint propagation, every t2 slice is coil''*rho');

% Check random SO(3) walks start at the identity orientation and remain finite
rng(1,'twister');
local_ensure_pool();
eulers=rwalk(5,1,1e-6);
result=test_true(result,'rwalk trajectory shape',isequal(size(eulers),[5 3])&&...
                 all(isfinite(eulers(:))),...
                 'rwalk() must return one finite Euler-angle triplet per trajectory point');
result=test_close(result,'rwalk initial orientation',euler2dcm(eulers(1,:)),eye(3),1e-14,1e-14,...
                  'the first random-walk orientation must reconstruct the starting reference frame');

end


function spin_system=local_liouvillian_system(dim)

% Create a quiet spherical-tensor Liouville descriptor
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.bas.formalism='sphten-liouv';
spin_system.bas.basis=zeros(dim,1);
spin_system.tols.liouv_zero=1e-12;
spin_system.tols.prop_chop=1e-12;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.small_matrix=10;

end


function spin_system=local_dipolar_system()

% Create a quiet two-spin descriptor with one close internuclear vector
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.comp.nspins=2;
spin_system.chem.parts={1:2};
spin_system.inter.coordinates={[0 0 0],[1 0 0]};
spin_system.inter.pbc={};
spin_system.inter.gammas=[2 3];
spin_system.inter.coupling.matrix=cell(2,2);
spin_system.tols.prox_cutoff=2;
spin_system.tols.dd_ncells=0;
spin_system.tols.hbar=1;
spin_system.tols.mu0=4*pi;

end


function [sys,inter]=local_isoswap_inputs()

% Create two-spin interaction structures with transferable and wiped terms
sys.isotopes={'1H','13C'};
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,1}=eye(3);
inter.coupling.matrix{1,2}=eye(3);
inter.coupling.matrix{2,1}=2*eye(3);
inter.coupling.matrix{2,2}=[];

end


function local_ensure_pool()

% Start a one-worker process pool for compact parfor utilities
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

end
