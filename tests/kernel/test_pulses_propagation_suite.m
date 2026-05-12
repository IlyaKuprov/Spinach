% Tests pulse-coordinate and propagation helpers. Syntax:
%
%                    result=test_pulses_propagation_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks RF coordinate Hessian round-trips, Iserles generators,
% Lie-step methods on a constant generator, and R-sequence phase/compiler
% invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_pulses_propagation_suite()

% State the propagation-helper target of the test
result=new_test_result('kernel/pulses_propagation_suite',...
                       'Pulse coordinate and propagation helpers',...
                       'pulse propagation helpers must preserve analytic coordinate and propagator identities.');

% Choose non-zero amplitudes away from the polar singularity
r=[1.2 2.3 3.4];
p=[0.2 -0.7 1.1];

% Use f=sum(r.^2)=sum(x.^2+y.^2), whose Cartesian Hessian is exactly 2I
Dr=2*r;
Dp=zeros(size(r));
Drr=2*eye(numel(r));
Drp=zeros(numel(r));
Dpr=zeros(numel(r));
Dpp=zeros(numel(r));

% Convert polar coordinates, gradients, and Hessians to Cartesian and back
[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp);
[r_back,p_back,Dr_back,Dp_back,Drr_back,Drp_back,Dpr_back,Dpp_back]=...
    cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy);
result=test_close(result,'polar-cartesian amplitude Hessian round-trip',r_back,r,1e-12,1e-12,...
                  'non-zero RF amplitudes survive polar-Cartesian coordinate round-trip');
result=test_close(result,'polar-cartesian phase Hessian round-trip',p_back,p,1e-12,1e-12,...
                  'RF phases away from branch cuts survive polar-Cartesian coordinate round-trip');
result=test_close(result,'polar-cartesian X gradient',Dx,2*x,1e-12,1e-12,...
                  'the gradient of RF power with respect to X is 2X');
result=test_close(result,'polar-cartesian Y gradient',Dy,2*y,1e-12,1e-12,...
                  'the gradient of RF power with respect to Y is 2Y');
result=test_close(result,'polar-cartesian Dxx Hessian',Dxx,2*eye(numel(r)),1e-12,1e-12,...
                  'the X-X Hessian block of RF power is 2I');
result=test_close(result,'polar-cartesian Dyy Hessian',Dyy,2*eye(numel(r)),1e-12,1e-12,...
                  'the Y-Y Hessian block of RF power is 2I');
result=test_close(result,'polar-cartesian Dxy Hessian',Dxy,zeros(numel(r)),1e-12,1e-12,...
                  'the X-Y Hessian block of RF power is zero');
result=test_close(result,'polar-cartesian Dyx Hessian',Dyx,zeros(numel(r)),1e-12,1e-12,...
                  'the Y-X Hessian block of RF power is zero');
result=test_close(result,'polar-cartesian amplitude-gradient round-trip',Dr_back,Dr,1e-10,1e-10,...
                  'amplitude gradients transform back by the chain rule');
result=test_close(result,'polar-cartesian phase-gradient round-trip',Dp_back,Dp,1e-10,1e-10,...
                  'phase gradients transform back by the chain rule');
result=test_close(result,'polar-cartesian Drr round-trip',Drr_back,Drr,1e-9,1e-9,...
                  'amplitude-amplitude Hessian blocks transform back by the second-order chain rule');
result=test_close(result,'polar-cartesian Drp round-trip',Drp_back,Drp,1e-9,1e-9,...
                  'amplitude-phase Hessian blocks transform back by the second-order chain rule');
result=test_close(result,'polar-cartesian Dpr round-trip',Dpr_back,Dpr,1e-9,1e-9,...
                  'phase-amplitude Hessian blocks transform back by the second-order chain rule');
result=test_close(result,'polar-cartesian Dpp round-trip',Dpp_back,Dpp,1e-9,1e-9,...
                  'phase-phase Hessian blocks transform back by the second-order chain rule');

% Check Iserles second-order and fourth-order product quadrature formulae
HL=[0 1;1 0];
HM=[0 -1i;1i 0];
HR=[1 0;0 -1];
dt=0.125;
H_ref=(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL);
result=test_close(result,'isergen second order',isergen(HL,[],HR,dt),H_ref,1e-15,1e-15,...
                  'second-order Iserles quadrature uses the endpoint average and commutator correction');
H_ref=(HL+4*HM+HR)/6+(1i*dt/12)*(HL*HR-HR*HL);
result=test_close(result,'isergen fourth order',isergen(HL,HM,HR,dt),H_ref,1e-15,1e-15,...
                  'fourth-order Iserles quadrature uses Simpson weights and a smaller commutator correction');

% Build a one-proton Hilbert-space spin system for Lie-step checks
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
L=operator(spin_system,'Lx',1);
rho=state(spin_system,'Lz',1);
dt=0.25;
exact=step(spin_system,L,rho,dt);
Lfun=@(~,~)L;

% Constant generators should reduce all low-order Lie methods to the exact exponential step
methods={'PWCL','LG2','LG4','RKMK4','LG4A'};
for n=1:numel(methods)
    rho_obs=iserstep(spin_system,{Lfun,0,methods{n}},rho,dt);
    result=test_close(result,['iserstep constant ' methods{n}],rho_obs,exact,1e-12,1e-12,...
                      'a state-independent generator has no time-ordering correction beyond the exact step');
end

% Check R-sequence phase construction and nutation calibration
[phases,pulse_amp,pulse_dur]=rsequence(1,4,1,1,1000,'180_pulse','homo_double_quantum_nucycle');
base=[pi/4;-pi/4;pi/4;-pi/4];
phase_ref=[base;-base];
result=test_close(result,'rsequence phases',phases,phase_ref,1e-15,1e-15,...
                  'R4_1 phase alternation is +/-pi*nu/N and the nucycle appends the inverted block');
result=test_close(result,'rsequence pulse amplitude',pulse_amp,4000*pi,1e-12,1e-12,...
                  'a 180-degree pulse over one quarter rotor period has nutation pi/duration');
result=test_close(result,'rsequence pulse duration',pulse_dur,1/(4*1000),1e-15,1e-15,...
                  'the R element duration is n rotor periods divided by N blocks');

% Check R-sequence compiler indexing when the RF amplitude is zero
S=pauli(2);
Sx=S.x;
Sy=S.y;
[P,T]=rseq_compiler(spin_system,zeros(2),Sx,Sy,[0;pi;0],0,0.1,'180_pulse');
result=test_close(result,'rseq_compiler index map',T,[1;2;1],1e-15,1e-15,...
                  'unique phases are compiled once and reused by their index map');
for n=1:numel(P)
    result=test_close(result,['rseq_compiler zero-amplitude propagator ' int2str(n)],P{n},speye(2),1e-15,1e-15,...
                      'with zero RF amplitude and zero drift each compiled propagator is the identity');
end

end


