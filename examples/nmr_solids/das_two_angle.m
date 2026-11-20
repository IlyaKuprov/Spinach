% Sequential double-angle spinning of an axial spin-1/2 CSA tensor
% in the high-field secular approximation. The axis jumps about the
% laboratory Y axis, carrying the crystallite rigidly with the rotor.
% The jump is instantaneous; neither the density nor the receiver is
% rotated in the laboratory spin frame.
%
% talos@spindynamics.org

function das_two_angle()

% Define the signed secular chemical shift in rad/s
rate=6000;
iso=2*pi*300;
tensor=2*pi*1200*diag([-1/2 -1/2 1]);
tilts=deg2rad([40 75]);
t_jump=0.43/rate;
t_after=(0:80)/(32*rate);
grid=load('leb_3ang_rank_5.mat','alphas','betas','gammas','weights');
traces=complex(zeros(2,numel(t_after)));
fine=complex(zeros(1,numel(t_after)));

% Keep each crystallite and rotor phase continuous across the jump
for n=1:numel(grid.weights)

    % Orient the crystallite once in the common rotor frame
    crystal=euler2dcm(grid.alphas(n),grid.betas(n),grid.gammas(n));
    Q=crystal*tensor*crystal';
    first=delay_phase(Q,tilts(1),rate,0,t_jump,iso);
    same=delay_phase(Q,tilts(1),rate,t_jump,t_jump+t_after,iso);
    jump=delay_phase(Q,tilts(2),rate,t_jump,t_jump+t_after,iso);
    continuous=delay_phase(Q,tilts(1),rate,0,t_jump+t_after,iso);
    assert(max(abs(first+same-continuous))<1e-11);
    weight=grid.weights(n);
    traces=traces+weight*exp(-1i*[first+same; first+jump]);

    % Direct midpoint time ordering uses explicit Cartesian rotor poses
    ref=0;
    left=0;
    for k=1:numel(t_after)

        % The first interval ends exactly at the tilt jump
        right=t_jump+t_after(k);
        tilt=tilts(2);
        if k==1
            tilt=tilts(1);
        end
        dt=(right-left)/64;
        for s=1:64
            mid=left+(s-0.5)*dt;
            pose=euler2dcm(0,tilt,0)*euler2dcm(-2*pi*rate*mid,0,0);
            lab=pose*Q*pose';
            ref=ref+dt*(iso+lab(3,3));
        end
        fine(k)=fine(k)+weight*exp(-1i*ref);
        left=right;

    end

end

% Verify lab-frame phase continuity and compare with direct evolution
assert(abs(traces(1,1)-traces(2,1))<1e-12);
error_max=max(abs(traces(2,:)-fine));
assert(error_max<1e-3&&norm(traces(1,:)-traces(2,:))>1e-3);
fprintf('DAS jumped trace versus direct midpoint: %.6g\n',error_max);

% Plot the real part of the same density-matrix coherence in both runs
kfigure(); plot(1e3*t_after,real(traces).');
kxlabel('Time after axis jump, ms');
kylabel('Real coherence'); kgrid;
klegend({'Equal axes','Axis jump'});

end

% Integrate lab-frame zz projections of the rigid rotor tensor
function phase=delay_phase(Q,tilt,rate,start,stop,iso)

% Five equispaced rotor samples resolve the complete rank-two harmonics
node=2*pi*(0:4)/5;
lab_z=[-sin(tilt)*cos(node); -sin(tilt)*sin(node); cos(tilt)*ones(size(node))];
samples=sum(lab_z.*(Q*lab_z),1);
coeff=fftshift(fft(samples))/5;
omega=2*pi*rate*[-2 -1 1 2].';

% Keep the DC contribution separate from the oscillating harmonics
osc=(exp(1i*omega*stop)-exp(1i*omega*start))./(1i*omega);
phase=real((iso+coeff(3))*(stop-start)+coeff([1 2 4 5])*osc);

end


