% Analytic conditional CCE limits, Hahn modulation and complex FID. Syntax:
%
%                         result=test_cce_analytic()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_analytic()

% State independent targets rather than comparing two expansion routes
result=new_test_result('kernel/cce_analytic','CCE analytic limits',...
                       'Zero bath, isolated phase and one-spin Hahn modulation.');

% Build a spin-half pair with explicitly prescribed angular-frequency terms
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2); inter.coupling.matrix{1,2}=zeros(3);
spin_system=create(sys,inter);
omega=1700; longitudinal=700; transverse=430; phase=230;
spin_system.inter.zeeman.matrix{1}(3,3)=phase;
spin_system.inter.zeeman.matrix{2}(3,3)=omega;
spin_system.inter.coupling.matrix{1,2}(1,3)=transverse;
spin_system.inter.coupling.matrix{1,2}(3,3)=longitudinal;
p=parameters(); p.times=linspace(0,0.012,81);

% Compare the complex conditional FID against independent two-by-two expm
fid=cce(spin_system,p); z=diag([1 -1])/2; x=[0 1;1 0]/2;
ha=omega*z+(longitudinal*z+transverse*x)/2;
hb=omega*z-(longitudinal*z+transverse*x)/2;
reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    t=p.times(k);
    reference(k)=exp(-1i*phase*t)*trace(expm(-1i*ha*t)*expm(1i*hb*t))/2;
end
result=test_close(result,'complex FID',fid.value,reference,1e-11,1e-11,'Independent analytic or dense reference.');

% Evaluate the exact single-spin Hahn formula for two non-collinear fields
p.sequence=struct('type','ideal','target',1,'placement','fraction','at',0.5,...
                  'axis','x','angle',pi);
echo=cce(spin_system,p); a=[transverse/2 0 omega+longitudinal/2];
b=[-transverse/2 0 omega-longitudinal/2];
modulation=norm(cross(a,b))^2/(norm(a)^2*norm(b)^2);
reference=1-2*modulation*sin(norm(a)*p.times(:)/4).^2.*sin(norm(b)*p.times(:)/4).^2;
result=test_close(result,'one-spin Hahn formula',echo.value,reference,1e-11,1e-11,'Independent analytic or dense reference.');
result=test_true(result,'zero time is unity',echo.value(1)==1,'Independent analytic or dense reference.');

% Verify isolated-centre phase appears once in an empty bath
p.bath_spins=[]; p.order=0; p.sequence=[];
empty=cce(spin_system,p);
result=test_close(result,'empty bath phase',empty.value,exp(-1i*phase*p.times(:)),1e-11,1e-11,'Independent analytic or dense reference.');

% Verify a disconnected unpolarized bath changes no central dynamics
p.bath_spins=2; p.order=1; spin_system.inter.coupling.matrix{1,2}=zeros(3);
uncoupled=cce(spin_system,p);
result=test_close(result,'uncoupled bath',uncoupled.value,empty.value,1e-11,1e-11,'Independent analytic or dense reference.');
metrics.fid_error=max(abs(fid.value-reference_fid(ha,hb,phase,p.times)));
metrics.hahn_error=max(abs(echo.value-reference));

% Report the measured metrics
disp(metrics);

end

% Supply all physics choices explicitly
function p=parameters()
p.central_spins=1; p.bath_spins=2; p.method='cce'; p.expansion='standard';
p.quantity='coherence'; p.order=1; p.central.alpha=[1;0]; p.central.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm');
p.bath_state=struct('mode','mixed_product'); p.meanfield='none'; p.sequence=[];
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
end

% Independently evaluate the analytic branch trace for metric reporting
function values=reference_fid(ha,hb,phase,times)
values=zeros(numel(times),1);
for k=1:numel(times)
    values(k)=exp(-1i*phase*times(k))*trace(expm(-1i*ha*times(k))*expm(1i*hb*times(k)))/2;
end
end


