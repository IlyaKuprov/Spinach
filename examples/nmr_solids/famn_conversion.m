% Finite-pulse 17O triple-quantum to central-transition conversion.
% Compares alternating-phase trains tuned for signed transfer or for
% transfer per unit time, using a small design and hold-out ensemble.
% This is a discrete design demonstration, not a powder convergence test.
%
% talos@spindynamics.org

function famn_conversion()

% Set the field, MAS rate, and finite RF amplitude
sys.magnet=2*pi*800e6/spin('1H');
sys.isotopes={'17O'};
sys.parallel={'processes',1};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
parameters.rate=20e3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=8;
parameters.spins={'17O'};
parameters.offset=0;
parameters.rframes={{'17O',2}};
parameters.masframe='rotor';
rf_amp=2*pi*140e3;
dwell=0.05e-6;

% Build a distinct hold-out set by changing both coupling and orientation
cq_hz=[7.1 7.1 6.8 7.4]*1e6;
angles=[0.2 0.7 0.4; 1.2 1.4 2.3; 0.9 1.1 0.5; 2.5 0.8 1.5];
prefix=cell(1,numel(cq_hz));
for n=1:numel(cq_hz)

    % Sample the rotor Hamiltonian and construct cumulative pulse propagators
    inter.coupling.matrix={eeqq2nqi(cq_hz(n),0.6,5/2,[0 0 0])};
    spin_system=basis(create(sys,inter),bas);
    parameters.orientation=angles(n,:);
    stack=rotor_stack(spin_system,parameters,'labframe');
    stack=cellfun(@full,stack,'UniformOutput',false);
    lx=full(operator(spin_system,'Lx',1));
    rotor_ticks=round(1/(parameters.rate*dwell));
    drift=interpft(cat(3,stack{:}),rotor_ticks,3);
    prefix{n}=repmat(eye(6),1,1,361,2);
    for k=1:360

        % Retain both phases on the same time grid
        midpoint=(drift(:,:,k)+drift(:,:,k+1))/2;
        for phase=1:2
            prefix{n}(:,:,k+1,phase)=expm(-1i*dwell*(midpoint+(-1)^(phase-1)*rf_amp*lx))*...
                                      prefix{n}(:,:,k,phase);
        end

    end

end

% Use matched phase conventions for the input and signed CT measurement
rho=sparse(5,2,1,6,6);
coil=centrans(6,'-');
nu_hz=-sys.magnet*spin('17O')/(2*pi);
fprintf('CT second-order isotropic correction: %.4f Hz\n',...
        quad_iso_shift(5/2,cq_hz(1),0.6,nu_hz,0.5,-0.5));
fprintf('Objective      Pulses  Duration/us   Train signed   Hold-out signed\n');

% Optimise one pulse, then append the opposite phase and retune all durations
for penalty=[0 1]

    % Optimise over the physical range 0.2 to 3 microseconds per pulse
    seed=20*ones(1,6);
    for n=1:6
        ticks=seed(1:n);

        % Cycle over pulse durations until no discrete change improves the score
        previous=[];
        while ~isequal(previous,ticks)
            previous=ticks;
            for k=1:n

                % Pick the best duration for this pulse without moving on ties
                trial=ticks;
                best=conversion(ticks,prefix(1:2),rho,coil,penalty);
                for candidate=4:60
                    trial(k)=candidate;
                    score=conversion(trial,prefix(1:2),rho,coil,penalty);
                    if score>best
                        ticks(k)=candidate;
                        best=score;
                    end
                end

            end
        end

        % Report the signed transfer independently of either design objective
        [~,train]=conversion(ticks,prefix(1:2),rho,coil,0);
        [~,held_out]=conversion(ticks,prefix(3:4),rho,coil,0);
        names={'Signed','Rate'};
        fprintf('%-12s %6d %12.3f %14.6f %17.6f\n',...
                names{penalty+1},n,1e6*dwell*sum(ticks),train,held_out);
        seed=[ticks repmat(ticks(end),1,6-n)];

    end

end

end

% Evaluate signed transfer and the selected optimisation objective
function [score,signal]=conversion(ticks,prefix,rho,coil,penalty)

% Propagate one matrix unit through the alternating pulse train
limits=[0 cumsum(ticks)]+1;
values=zeros(size(prefix));
for n=1:numel(prefix)

    % Compose finite pulse propagators in time order
    P=eye(6);
    for k=1:numel(ticks)
        phase=1+mod(k-1,2);
        P=prefix{n}(:,:,limits(k+1),phase)*prefix{n}(:,:,limits(k),phase)'*P;
    end
    values(n)=trace(coil'*P*rho*P');

end

% Sign-sensitive maximum or squared magnitude normalised by train length squared
overlap=mean(values);
signal=real(overlap);
if penalty==0
    score=signal;
else
    score=abs(overlap/sum(ticks))^2;
end

end


