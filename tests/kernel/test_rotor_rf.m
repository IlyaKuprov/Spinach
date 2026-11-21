% Test: kernel/rotor_rf | Rotor RF map
% Tests the RF gain map, its Jacobian, and linked catalog rows.
% Syntax: result=test_rotor_rf()
%
% Talos <talos@spindynamics.org>

function result=test_rotor_rf()

% Apply one complex sample at a known column-major address
result=new_test_result('kernel/rotor_rf','Rotor RF map','RF gains act on local quadrature pairs.');
wave=zeros(4,3);
wave(3,2)=1;
gain=repmat([1-1i;2+3i],1,3);
[out,J]=rotor_rf(wave,gain);
expected=zeros(size(wave));
expected(3:4,2)=[2;3];
result=test_close(result,'map',out,expected,0,0,'The complex multiplier acts on one pair.');
result=test_close(result,'Jacobian column',J(:,7),expected(:),0,0,'No time or channel mixing.');
result=test_true(result,'sparse real Jacobian',issparse(J)&&isreal(J),'The Jacobian is sparse and real.');

% Combine electronics and a field map, leaving one input quadrature frozen
wave=reshape(1:12,4,3)/10;
direction=cos(wave);
direction(1,2)=0;
[filtered,F]=firf(wave,[0.8 0.4]);
[~,J]=rotor_rf(filtered,gain);
delta=1e-5;
plus=rotor_rf(firf(wave+delta*direction,[0.8 0.4]),gain);
minus=rotor_rf(firf(wave-delta*direction,[0.8 0.4]),gain);
result=test_close(result,'composed derivative',J*F*direction(:),...
                  (plus(:)-minus(:))/(2*delta),1e-9,0,'The map follows the electronics Jacobian.');

% Compare linked and Cartesian averages with distinguishable drift-map pairs
control.rho_init={1};
control.ndrifts=2;
control.pwr_levels=1;
control.offsets={};
control.phase_cycle=zeros(1,0);
control.distortion={@no_dist;@no_dist};
control.budget=Inf;
control.ens_corrs={'dist_drift'};
control.drift_weights=[1;3];
[linked,~,weights]=ens_catalog(control);
control.ens_corrs={};
[independent,~,cart_weights]=ens_catalog(control);
response=[1 3; 5 9];
linked_val=weights.'*response(sub2ind([2 2],linked(:,2),linked(:,6)));
cart_val=cart_weights.'*response(sub2ind([2 2],independent(:,2),independent(:,6)));
result=test_close(result,'linked average',linked_val,7,0,0,'Only the matched pair contributes.');
result=test_close(result,'Cartesian average',cart_val,5.75,0,0,'Independent maps mix different drift cases.');

% Exercise the linked maps through the public objective with a frozen input
sys.isotopes={'1H'};
sys.magnet=1;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=test_spin_system(sys,struct(),bas);
S=pauli(2);
control=struct();
control.isotopes={'1H'};
control.channels=[1;1];
control.operators={S.x,S.y};
control.rho_init={S.z};
control.rho_targ={S.x};
control.drifts={{S.z},{-2*S.z}};
control.pulse_dt=[0.1 0.17];
control.pwr_levels=1;
control.penalties={'none'};
control.p_weights=0;
control.ens_corrs={'dist_drift'};
control.drift_weights=[1 3];
control.freeze=logical([1 0;0 0]);
control.distortion={@(w)firf(w,[0.8 0.4]),@(w)rotor_rf(w,[1+1i 0.5]);...
                    @(w)firf(w,[0.8 0.4]),@(w)rotor_rf(w,[0.2-1i 1.1])};
local_system=optimcon(spin_system,control);
wave=[2 -1;1 3];
[~,fid,grad]=grape_xy(wave,local_system);

% Hand-average independent one-case objective evaluations
expected=0;
for n=1:2

    % Retain only one drift and its matching electronics and map
    member=rmfield(control,{'ens_corrs','drift_weights'});
    member.drifts=control.drifts(n);
    member.distortion=control.distortion(n,:);
    [~,val]=grape_xy(wave,optimcon(spin_system,member));
    expected=expected+control.drift_weights(n)*val(1)/4;

end
result=test_close(result,'linked physical objective',fid(1),expected,1e-12,0,...
                  'The public objective uses the matched subset of real trajectories.');
unlinked=rmfield(control,'ens_corrs');
[~,mixed]=grape_xy(wave,optimcon(spin_system,unlinked));
result=test_true(result,'unlinked physical objective',abs(fid(1)-mixed(1))>1e-3,...
                 'Deliberately mismatched drift-map combinations change the observable.');

% Compare the complete phase-space pullback on unfrozen coordinates
direction=[0 0.7;-0.3 0.2];
[~,plus]=grape_xy(wave+delta*direction,local_system);
[~,minus]=grape_xy(wave-delta*direction,local_system);
result=test_close(result,'linked objective derivative',sum(grad(:,:,1).*direction,'all'),...
                  (plus(1)-minus(1))/(2*delta),1e-9,0,'Finite differences include both distortion stages.');
result=test_close(result,'frozen quadrature',grad(1,1,1),0,0,0,'Input freezing survives the pullback.');

end


