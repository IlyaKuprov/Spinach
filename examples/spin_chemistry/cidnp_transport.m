% Nuclear polarisation transported from a radical pair to a diamagnetic
% singlet-channel product. The electrons are traced out by atom matching;
% the triplet channel leads to a tracked spin-free sink. The effective
% electron offsets and hyperfine coupling below are in rad/s, and the
% recombination rates are in inverse seconds. This is a small illustrative
% model, not a fitted molecular system.
%
% ilya.kuprov@weizmann.ac.il

function cidnp_transport()

% Radical pair, nuclear product, and spin-free sink
sys.magnet=1; sys.isotopes={'E','E','1H','1H'};
inter.chem.parts={1:3,4,[]}; inter.chem.concs=[1 0 0];
inter.chem.reactions={...
    struct('reactants',1,'products',2,'matching',[3 4],'rate',2,...
           'selector',{{'singlet',[1 2]}}),...
    struct('reactants',1,'products',3,'matching',zeros(0,2),'rate',3,...
           'selector',{{'triplet',[1 2]}})};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Effective electron offsets and isotropic hyperfine interaction
H=2*operator(spin_system,'Lz',1)+1.7*operator(spin_system,'Lz',2)...
  +1.2*operator(spin_system,{'Lz','Lz'},{1,3})...
  +0.6*operator(spin_system,{'L+','L-'},{1,3})...
  +0.6*operator(spin_system,{'L-','L+'},{1,3});
K=kinetics(spin_system,'report');

% Unit-concentration electron singlet with an unpolarised nucleus
rho=4*singlet(spin_system,1,2);
coil=coil_state(spin_system,'Lz',4,'exact');
time_grid=linspace(0,3,301); dt=time_grid(2)-time_grid(1);
traj=evolution(spin_system,H+1i*K,[],rho,dt,numel(time_grid)-1,'trajectory');

% Channel populations and transported nuclear polarisation
concs=zeros(numel(time_grid),3);
for n=1:numel(time_grid)
    concs(n,:)=real(chem_concs(spin_system,traj(:,n)));
end
signal=real(coil'*traj);
fprintf('maximum concentration-balance error %.6g\n',max(abs(sum(concs,2)-1)));
fprintf('final product nuclear polarisation %.6g\n',signal(end));

% Show nuclear arrival without the unobserved electron coordinates
kfigure(); subplot(2,1,1); plot(time_grid,concs); kgrid;
kxlabel('time, seconds'); kylabel('concentration');
klegend({'radical pair','singlet product','triplet sink'});
subplot(2,1,2); plot(time_grid,signal); kgrid;
kxlabel('time, seconds'); kylabel('product nuclear polarisation');

end


