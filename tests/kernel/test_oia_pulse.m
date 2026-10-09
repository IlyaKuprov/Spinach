% Tests offset-independent adiabaticity pulses. Syntax:
%
%                    result=test_oia_pulse()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that oia_pulse() reproduces chirp_pulse() for a
% constant amplitude function, that its numerically integrated fre-
% quency sweeps match the closed-form sweeps in Table 1 of Tannus
% and Garwood (JMR A 120, 133 (1996)), and that every amplitude
% function from that table, as well as an arbitrary smooth one,
% inverts a single spin to better than 98% at 17 offsets across
% the central 80% of the sweep bandwidth and leaves it at least
% 90% longitudinal at offsets 30% beyond the sweep edges; the
% transition bands are not tested.
%
% ilya.kuprov@weizmann.ac.il

function result=test_oia_pulse()

% Announce the test target
fprintf('TESTING: Offset-independent adiabaticity pulses\n');

% State the physics target of the test
result=new_test_result('kernel/oia_pulse',...
                       'Offset-independent adiabaticity pulses',...
                       'the frequency sweep that follows the squared amplitude must invert uniformly inside the bandwidth.');

% Pulse parameters from the paper: 50 kHz sweep in 2 ms
npts=1000; dur=2e-3; bwidth=50e3;

% Constant amplitude must reproduce the square chirp
[Cx_oia,Cy_oia]=oia_pulse(npts,dur,bwidth,@(tau)ones(size(tau)));
[Cx_chirp,Cy_chirp]=chirp_pulse(npts,dur,bwidth,0,'smoothed');
result=test_close(result,'square chirp coincidence',[Cx_oia; Cy_oia],[Cx_chirp; Cy_chirp],1e-9,1e-9,...
                  'a constant amplitude function must give the square envelope chirp of chirp_pulse()');

% Closed-form sweeps from Table 1 of the paper, 1% edge truncation
beta_hs=asech(0.01); beta_gs=sqrt(2*log(100)); beta_lz=99;
am_funs={@(tau)sech(beta_hs*tau),...
         @(tau)exp(-beta_gs^2*tau.^2/2),...
         @(tau)1./(1+beta_lz*tau.^2),...
         @(tau)(1+cos(pi*tau))/2};
fm_funs={@(tau)tanh(beta_hs*tau)/tanh(beta_hs),...
         @(tau)erf(beta_gs*tau)/erf(beta_gs),...
         @(tau)(tau./(1+beta_lz*tau.^2)+atan(sqrt(beta_lz)*tau)/sqrt(beta_lz))/...
               (1/(1+beta_lz)+atan(sqrt(beta_lz))/sqrt(beta_lz)),...
         @(tau)tau+(4/(3*pi))*sin(pi*tau).*(1+cos(pi*tau)/4)};
fm_names={'HS','Gauss','Lorentz','Hanning'};
tau=linspace(-1,1,20001);
for n=1:numel(am_funs)
    [~,~,~,~,~,~,frqs]=oia_pulse(numel(tau),dur,bwidth,am_funs{n});
    result=test_close(result,[fm_names{n} ' sweep closed form'],2*frqs/bwidth,fm_funs{n}(tau),1e-6,1e-6,...
                      'the integral of the squared amplitude must match the analytical sweep in Table 1');
end

% Single proton in Liouville space
sys.magnet=14.1; sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);

% Operators and the initial state
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');
rho0=state(spin_system,'Lz','1H');
rho0=rho0/norm(rho0,2);

% Amplitude functions from Table 1 and an arbitrary smooth one
am_funs={@(tau)1./(1+99*tau.^2),...
         @(tau)sech(asech(0.01)*tau),...
         @(tau)exp(-log(100)*tau.^2),...
         @(tau)(1+cos(pi*tau))/2,...
         @(tau)sech(asech(0.01)*tau.^8),...
         @(tau)1-abs(sin(pi*tau/2)).^40,...
         @(tau)(1-tau.^2).*(1.3+sin(2*tau))};
am_names={'Lorentz','HS','Gauss','Hanning','HS8','Sin40','hand drawn'};

% Offsets inside and outside the sweep bandwidth
offs_in=linspace(-0.8,0.8,17)*bwidth/2;
offs_out=[-1.3 1.3]*bwidth/2;
for n=1:numel(am_funs)

    % Make the pulse
    [Cx,Cy,durs]=oia_pulse(npts,dur,bwidth,am_funs{n});

    % Inversion profile inside the bandwidth
    mz_in=zeros(size(offs_in));
    for k=1:numel(offs_in)
        rho=shaped_pulse_xy(spin_system,2*pi*offs_in(k)*Lz,{Lx,Ly},{Cx,Cy},durs,rho0,'expv-pwc');
        mz_in(k)=real(rho0'*rho);
    end
    result=test_true(result,[am_names{n} ' inversion inside the bandwidth'],all(mz_in<-0.98),...
                     'every offset inside 80% of the sweep must be inverted to better than 98%');

    % Magnetisation outside the bandwidth
    mz_out=zeros(size(offs_out));
    for k=1:numel(offs_out)
        rho=shaped_pulse_xy(spin_system,2*pi*offs_out(k)*Lz,{Lx,Ly},{Cx,Cy},durs,rho0,'expv-pwc');
        mz_out(k)=real(rho0'*rho);
    end
    result=test_true(result,[am_names{n} ' no inversion outside the bandwidth'],all(mz_out>0.9),...
                     'offsets 30% beyond the sweep edge must retain at least 90% of the longitudinal magnetisation');

end

end


