% A+B to C with additive and product spin-arrival closures. Each reactant
% carries longitudinal polarisation; only product closure creates their
% cross-reactant two-spin order. RKMK4 is compared with ode45 at relative
% and absolute tolerances 1e-12 and 1e-14, respectively.
%
% ilya.kuprov@weizmann.ac.il

function bimolecular_closures()

% Two one-proton reactants and a two-proton product
sys.magnet=1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1,2,3:4}; inter.chem.concs=[0.7 0.3 0];
inter.chem.reactions={struct('reactants',[1 2],'products',3,...
                            'matching',[1 3;2 4],'rate',50)};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Concentration-weighted initial state and unweighted product detector
rho0=unit_state(spin_system)+0.1*state(spin_system,'Lz',1)...
                            +0.2*state(spin_system,'Lz',2);
coil=coil_state(spin_system,{'Lz','Lz'},{3,4},'exact');
time_grid=linspace(0,0.02,41); dt=time_grid(2)-time_grid(1);
closures={'additive','product'}; signal=zeros(2,numel(time_grid));

% Compare the two nonlinear closures with independently stepped references
for n=1:2
    spin_system.chem.reactions{1}.closure=closures{n};
    K=kinetics(spin_system,'report'); rho=rho0;
    signal(n,1)=real(coil'*rho);
    for k=2:numel(time_grid)
        rho=step(spin_system,{@(t,y)1i*K(t,y),time_grid(k-1),'RKMK4'},rho,dt);
        signal(n,k)=real(coil'*rho);
    end
    [~,reference]=ode45(@(t,y)K(t,y)*y,[time_grid(1) time_grid(end)],...
                        full(rho0),odeset('RelTol',1e-12,'AbsTol',1e-14));
    fprintf('%s closure: endpoint error against ode45 %.6g\n',...
            closures{n},norm(rho-reference(end,:)'));
end

% Product two-spin order is absent only in the additive closure
kfigure(); plot(time_grid,signal'); kgrid;
kxlabel('time, seconds'); kylabel('product two-spin order');
klegend(closures);

end


