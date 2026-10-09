% Reversible capture A+Q to S, where A carries a proton, Q is a spin-free
% quencher, and S is a spin-free product. Reverse reaction creates an
% unpolarised proton. All populations are unit coordinates in one state;
% A+Q+2*S is the conserved atom-equivalent concentration.
%
% ilya.kuprov@weizmann.ac.il

function spinless_sink_network()

% Reaction network and initially empty sink
sys.magnet=1; sys.isotopes={'1H'};
inter.chem.parts={1,[],[]}; inter.chem.concs=[0.7 0.3 0];
inter.chem.reactions={...
    struct('reactants',[1 2],'products',3,'matching',zeros(0,2),'rate',5),...
    struct('reactants',3,'products',[1 2],'matching',zeros(0,2),'rate',1)};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial polarisation and concentration-weighted populations
rho=unit_state(spin_system)+0.2*state(spin_system,'Lz',1);
coil=coil_state(spin_system,'Lz',1,'exact');
K=kinetics(spin_system,'report');
time_grid=linspace(0,1,201); dt=time_grid(2)-time_grid(1);
concs=zeros(numel(time_grid),3); signal=zeros(size(time_grid));
concs(1,:)=chem_concs(spin_system,rho); signal(1)=real(coil'*rho);

% Propagate concentration and spin order together
for n=2:numel(time_grid)
    rho=step(spin_system,{@(t,y)1i*K(t,y),time_grid(n-1),'RKMK4'},rho,dt);
    concs(n,:)=real(chem_concs(spin_system,rho)); signal(n)=real(coil'*rho);
end
fprintf('maximum atom-equivalent balance error %.6g\n',max(abs(concs*[1;1;2]-1)));

% Compare species populations with the surviving magnetic signal
kfigure(); subplot(2,1,1); plot(time_grid,concs); kgrid;
kxlabel('time, seconds'); kylabel('concentration'); klegend({'A','Q','S'});
subplot(2,1,2); plot(time_grid,signal); kgrid;
kxlabel('time, seconds'); kylabel('proton polarisation');

end


