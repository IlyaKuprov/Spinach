% Accuracy test for Lie-group product quadratures as a function of
% discretisation step in the E1000B Veshtort-Griffin pulse.
%
% Calculation time: seconds
%
% a.graham@soton.ac.uk
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function product_quadratures_1()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Zeeman interactions
inter.zeeman.scalar=num2cell(linspace(-4,4,31));

% Couplings
inter.coupling.scalar=cell(31);
for n=1:30
    inter.coupling.scalar{n,n+1}=10;
end

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'nmr');

% Hamiltonian superoperator
H=hamiltonian(spin_system);

% Control operator
Lx=operator(spin_system,'Lx','1H');

% Pulse duration
duration=1e-2;

% Accurate reference calculation
rho_ref=state(spin_system,'Lz','1H');
np=2000; dt=duration/(np-1);
amps=vg_pulse('E1000B',2*np-1,duration);
for n=1:2:(numel(amps)-2)
    rho_ref=step(spin_system,{H+amps(n+0)*Lx,...
                              H+amps(n+1)*Lx,...
                              H+amps(n+2)*Lx},rho_ref,dt);
end

% Plotting
kfigure(); scale_figure([1.50 0.60]);
time_axis=linspace(0,duration,2*np-1)';
subplot(1,2,1); plot(time_axis,amps/(2*pi)); 
xlim tight; ylim padded; kgrid;
kylabel('nutation freq., Hz');
kxlabel('time, seconds');

% Benchmark arrays
npi=100:100:1000; bench=zeros(numel(npi),4);

% Benchmarking loop
parfor k=1:numel(npi)

    % Left point integrator
    rho_left=state(spin_system,'Lz','1H');
    np=npi(k); dt=duration/(np-1);
    amps=vg_pulse('E1000B',np,duration);
    for n=1:(np-1)
        rho_left=step(spin_system,H+amps(n)*Lx,rho_left,dt);
    end

    % Midpoint integrator
    rho_mid=state(spin_system,'Lz','1H');
    np=npi(k); dt=duration/(np-1);
    amps=vg_pulse('E1000B',np,duration);
    for n=1:(np-1)
        rho_mid=step(spin_system,H+0.5*(amps(n+0)+...
                                        amps(n+1))*Lx,rho_mid,dt);
    end

    % Two-point integrator
    rho_two=state(spin_system,'Lz','1H');
    np=npi(k); dt=duration/(np-1);
    amps=vg_pulse('E1000B',np,duration);
    for n=1:(np-1)
        rho_two=step(spin_system,{H+amps(n+0)*Lx,...
                                  H+amps(n+1)*Lx},rho_two,dt);
    end

    % Three-point integrator
    rho_thr=state(spin_system,'Lz','1H');
    np=npi(k); dt=duration/(np-1);
    amps=vg_pulse('E1000B',2*np-1,duration);
    for n=1:2:(numel(amps)-2)
        rho_thr=step(spin_system,{H+amps(n+0)*Lx,...
                                  H+amps(n+1)*Lx,...
                                  H+amps(n+2)*Lx},rho_thr,dt);
    end
    
    % Benchmarks
    bench(k,:)=[norm(rho_ref-rho_left)/norm(rho_ref) ...
                norm(rho_ref-rho_mid)/norm(rho_ref) ...
                norm(rho_ref-rho_two)/norm(rho_ref) ...
                norm(rho_ref-rho_thr)/norm(rho_ref)];
    
end

% Plotting
subplot(1,2,2);
plot(npi',bench); axis tight;
ylim([1e-9 0.2]);
set(gca,'YScale','log'); 
set(gca,'XScale','log'); kgrid;
set(gca,'xminorgrid','off');
set(gca,'yminorgrid','off');
klegend({'LP','MP','LG-2','LG-4'},...
         'Location','southwest');
kxlabel('number of points in the time grid');
kylabel('$\|$difference$\|$/$\|$exact$\|$');
set(gca,'XTick',[100 200 400 800]);
set(gca,'YTick',10.^(-9:2:-1));

end

