% Tests Shinnar-Le Roux selective excitation pulse design. Syntax:
%
%                         result=test_slr_pulse()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks waveform units and shape, independent two-level
% propagation, excitation profile selectivity, production-path shaped
% pulse propagation, and representative input validation failures.
%
% ilya.kuprov@weizmann.ac.il

function result=test_slr_pulse()

% Announce the test target
fprintf('TESTING: Shinnar-Le Roux pulse design\n');

% State the selective pulse target of the test
result=new_test_result('kernel/slr_pulse',...
                       'Shinnar-Le Roux pulse design',...
                       'an SLR waveform must produce the requested flip and a selective excitation profile.');

% Define a representative selective excitation design
npts=64;
dur=4e-3;
tbw=4;
flip_angle=pi/2;
pass_rip=0.01;
stop_rip=0.01;

% Generate the production waveform
[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip);

% Check the output dimensions and finiteness
result=test_true(result,'SLR output dimensions',...
                 isequal(size(Cx),[1 npts])&&isequal(size(Cy),[1 npts])&&...
                 isequal(size(durs),[1 npts])&&isequal(size(amps),[1 npts])&&...
                 isequal(size(phis),[1 npts]),...
                 'all waveform outputs must be 1 x npts row vectors');
result=test_true(result,'SLR finite outputs',...
                 all(isfinite([Cx Cy durs amps phis])),...
                 'all generated waveform coordinates must be finite');
result=test_close(result,'SLR duration sum',sum(durs),dur,1e-15,1e-14,...
                  'uniform pulse slice durations must sum to the requested total duration');

% Check Cartesian and polar coordinate consistency
result=test_close(result,'SLR polar X consistency',amps.*cos(phis),Cx,1e-10,1e-14,...
                  'polar amplitude and phase must reproduce the Cartesian X control');
result=test_close(result,'SLR polar Y consistency',amps.*sin(phis),Cy,1e-10,1e-14,...
                  'polar amplitude and phase must reproduce the Cartesian Y control');

% Build independent spin-half matrices
Lx=[0 1;1 0]/2;
Ly=[0 -1i;1i 0]/2;
Lz=[1 0;0 -1]/2;

% Propagate the waveform directly at zero offset
rho=Lz;
for slice=1:npts
    U=expm(-1i*(Cx(slice)*Lx+Cy(slice)*Ly)*durs(slice));
    rho=U*rho*U';
end

% Check the requested flip and Spinach phase convention
result=test_close(result,'SLR direct centre flip',rho,-Ly,1e-11,1e-11,...
                  'a positive-X pi/2 control under exp(-1i*H*t) rotates Lz to -Ly');

% Check an interior value of the documented flip-angle interval
check_flip=pi/6;
[small_x,small_y,small_durs]=slr_pulse(npts,dur,tbw,check_flip,pass_rip,stop_rip);
rho=Lz;
for slice=1:npts
    U=expm(-1i*(small_x(slice)*Lx+small_y(slice)*Ly)*small_durs(slice));
    rho=U*rho*U';
end
flip_ref=cos(check_flip)*Lz-sin(check_flip)*Ly;
result=test_close(result,['SLR centre flip ' num2str(check_flip/pi) ' pi'],...
                  rho,flip_ref,1e-11,1e-11,...
                  'the beta-polynomial target must produce the documented interior flip');

% Build a dense independent Cayley-Klein frequency sweep
freq_grid=linspace(-0.5,0.5,16385);
[transverse,longitudinal,unit_err]=ck_profile(Cx,Cy,durs,freq_grid);

% Reconstruct the documented transition edges independently
beta_pass=sqrt(pass_rip/2);
beta_stop=stop_rip/sqrt(2);
log_pass=log10(beta_pass);
log_stop=log10(beta_stop);
trans_measure=(5.309e-3*log_pass^2+7.114e-2*log_pass-4.761e-1)*log_stop+...
              (-2.66e-3*log_pass^2-5.941e-1*log_pass-4.278e-1);
pass_edge=(tbw-trans_measure)/(2*npts);
stop_edge=(tbw+trans_measure)/(2*npts);
pass_mask=(abs(freq_grid)<=0.8*pass_edge);
stop_mask=(abs(freq_grid)>=1.2*stop_edge);

% Check unitarity and the selective excitation profile
result=test_true(result,'SLR Cayley-Klein unitarity',unit_err<2e-12,...
                 'every frequency response must remain unitary throughout propagation');
result=test_true(result,'SLR passband excitation',...
                 min(transverse(pass_mask))>0.99&&max(abs(longitudinal(pass_mask)))<0.10,...
                 'the interior passband must produce near-complete transverse excitation');
result=test_true(result,'SLR stopband suppression',...
                 max(transverse(stop_mask))<0.03&&max(abs(longitudinal(stop_mask)-1))<1e-3,...
                 'the exterior stopband must retain longitudinal magnetisation');

% Check that a smaller flip is designed before the nonlinear inverse transform
[small_trans,small_long,small_unit]=ck_profile(small_x,small_y,small_durs,freq_grid);
result=test_true(result,'SLR small-flip unitarity',small_unit<2e-12,...
                 'the smaller-flip frequency response must remain unitary');
result=test_true(result,'SLR small-flip passband',...
                 max(abs(small_trans(pass_mask)-sin(check_flip)))<0.025&&...
                 max(abs(small_long(pass_mask)-cos(check_flip)))<0.015,...
                 'the smaller flip must retain its requested selective passband response');
result=test_true(result,'SLR small-flip stopband',...
                 max(small_trans(stop_mask))<0.01&&...
                 max(abs(small_long(stop_mask)-1))<5e-5,...
                 'the smaller flip must retain stopband suppression');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Apply the generated controls through the production shaped-pulse path
Lx_op=operator(spin_system,'Lx',1);
Ly_op=operator(spin_system,'Ly',1);
Lz_state=state(spin_system,'Lz',1);
Ly_state=state(spin_system,'Ly',1);
rho_obs=shaped_pulse_xy(spin_system,0*Lx_op,{Lx_op,Ly_op},...
                        {Cx,Cy},durs,Lz_state,'expm-pwc');

% Check production-path sign and units
result=test_close(result,'SLR shaped-pulse centre flip',rho_obs,-Ly_state,1e-10,1e-10,...
                  'the generated rad/s controls and second durations must produce the requested flip');

% Check representative consistency failures
result=test_true(result,'SLR odd sample count rejected',...
                 throws_with(@()slr_pulse(63,dur,tbw,flip_angle,pass_rip,stop_rip),'even'),...
                 'the linear-phase design requires an even sample count');
result=test_true(result,'SLR nonpositive duration rejected',...
                 throws_with(@()slr_pulse(npts,0,tbw,flip_angle,pass_rip,stop_rip),'positive'),...
                 'pulse duration must be strictly positive');
result=test_true(result,'SLR invalid ripple rejected',...
                 throws_with(@()slr_pulse(npts,dur,tbw,flip_angle,1,stop_rip),'strictly between'),...
                 'ripple targets must lie strictly between zero and one');
result=test_true(result,'SLR excessive flip rejected',...
                 throws_with(@()slr_pulse(npts,dur,tbw,pi,pass_rip,stop_rip),'pi/2'),...
                 'selective excitation flip angles must not exceed pi/2');
result=test_true(result,'SLR infeasible transition rejected',...
                 throws_with(@()slr_pulse(8,dur,0.1,flip_angle,pass_rip,stop_rip),'transition band'),...
                 'sampling and bandwidth must accommodate the requested transition');

end

% Check that a function call throws an expected error
function verdict=throws_with(fun_handle,message_text)
try
    fun_handle();
    verdict=false;
catch err
    verdict=contains(err.message,message_text);
end
end

% Evaluate a waveform with an independent Cayley-Klein propagator
function [transverse,longitudinal,unit_err]=ck_profile(Cx,Cy,durs,freq_grid)
offsets=2*pi*freq_grid/durs(1);
alpha=ones(size(freq_grid));
beta=zeros(size(freq_grid));
unit_err=0;
for slice=1:numel(Cx)
    rot_rate=sqrt(Cx(slice)^2+Cy(slice)^2+offsets.^2);
    sine_term=sin(rot_rate*durs(slice)/2)./rot_rate;
    sine_term(rot_rate==0)=durs(slice)/2;
    step_alpha=cos(rot_rate*durs(slice)/2)-1i*offsets.*sine_term;
    step_beta=-1i*(Cx(slice)+1i*Cy(slice))*sine_term;
    old_alpha=alpha;
    alpha=step_alpha.*alpha-conj(step_beta).*beta;
    beta=step_beta.*old_alpha+conj(step_alpha).*beta;
    unit_err=max(unit_err,max(abs(abs(alpha).^2+abs(beta).^2-1)));
end
transverse=abs(2*conj(alpha).*beta);
longitudinal=abs(alpha).^2-abs(beta).^2;
end

