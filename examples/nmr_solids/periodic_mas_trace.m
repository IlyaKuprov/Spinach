% Compares the periodic Hilbert-space MAS FID of one tilted
% crystallite with direct propagation on successively finer rotor
% grids. The acquisition dwell remains one third of a rotor period.
%
% talos@spindynamics.org

function periodic_mas_trace()

% One carbon with a non-axial chemical shift tensor
sys.magnet=9.4;
sys.isotopes={'13C'};
sys.parallel={'processes',1};
inter.zeeman.eigs={[-42 -14 62]};
inter.zeeman.euler={[0.3 0.7 -0.5]};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotor and acquisition settings
parameters.rate=6000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.spins=sys.isotopes;
parameters.offset=0;
parameters.rframes={};
parameters.grid='single_crystal';
parameters.orientation=[0 0 0];
parameters.masframe='rotor';
parameters.rho0=state(spin_system,'Lx','13C');
parameters.coil=coil_state(spin_system,'L+','13C');
parameters.tick_rf={};
parameters.start_tick=1;
parameters.dwell=1/(3*parameters.rate);
parameters.npoints=13;

% Keep the acquisition grid fixed while refining the rotor left endpoints
ranks=[4 13 40 121];
traces=complex(zeros(numel(ranks),parameters.npoints));
for n=1:numel(ranks)

    % Acquire with the periodic callback, then walk the independently built stack
    parameters.max_rank=ranks(n);
    traces(n,:)=singlerot(spin_system,@periodic_acquire,parameters,'nmr');
    drift=rotor_stack(spin_system,parameters,'nmr');
    tick_time=1/(parameters.rate*numel(drift));
    sample_stride=numel(drift)/3;
    rho=parameters.rho0;
    ref=complex(zeros(1,parameters.npoints));
    for k=1:numel(ref)
        ref(k)=trace(parameters.coil'*rho);
        for s=1:sample_stride
            idx=mod(-((k-1)*sample_stride+s-1),numel(drift))+1;
            U=expm(-1i*full(drift{idx})*tick_time);
            rho=U*rho*U';
        end
    end
    assert(norm(traces(n,:)-ref)<1e-10);

end

% Verify convergence separately from the equal-grid route comparison
errors=vecnorm(diff(traces,1,1),2,2);
assert(all(diff(errors)<0));
disp(errors);

% Plot the complex trace at the finest rotor resolution
time=1e3*parameters.dwell*(0:parameters.npoints-1);
kfigure(); plot(time,real(traces(end,:)),time,imag(traces(end,:)));
kxlabel('Acquisition time, ms');
kylabel('Receiver amplitude'); kgrid;
klegend({'Real','Imaginary'});

end


