% Proton-pair contribution to the Hahn echo of two electron spins:
% Fig. 3(b) of Chen et al. (https://doi.org/10.1063/5.0303171).
% The two polarised triplets select the branches in Eq. (29). In the
% antiparallel proton subspace the flip-flop matrix element is C and
% the diagonal energy difference is D+E or D-E. Common Zeeman and
% exchange energies do not enter this pair echo.
% Syntax:
%
%               [loss,axes_hz]=gcce_two_centres()
%
% loss is 1-|L| at total echo time 40 us, with one row per flip-flop
% coupling and one column per hyperfine detuning. The corresponding
% signed axes in Hz are returned in axes_hz.flip and axes_hz.detuning.
%
% ilya.kuprov@weizmann.ac.il

function [loss,axes_hz]=gcce_two_centres()

% Electron pair and antiparallel proton pair in the secular Hamiltonian
sys.isotopes={'E','E','1H','1H'};
sys.magnet=0.01;
sys.parallel={'processes',1};
sys.output='hush';
inter.zeeman.scalar={0,0,hz2ppm(5000,sys.magnet,'1H'),0};
inter.coupling.matrix=cell(4);

% Refocus the double-quantum coherence between the two polarised triplets
parameters.observed=[1 2];
parameters.transition=struct('spins',[1 2],'ms',[0.5 -0.5;0.5 -0.5]);
parameters.times=40e-6;
parameters.cce.order=2;
parameters.cce.approximation='generalized';
parameters.cce.basis.formalism='zeeman-hilb';
parameters.cce.bath_state=struct('mode','explicit_product',...
                                 'rho',{{diag([1 0]),diag([0 1])}});

% Sweep the independent quadrant of the pair-correlation surface
flip=0:500:4000;
detuning=0:2500:100000;
loss=zeros(numel(flip),numel(detuning));
for n=2:numel(flip)

    % Secular proton dipolar coupling with flip-flop matrix element C
    inter.coupling.matrix{3,4}=diag([2 2 -4])*flip(n);
    for k=2:numel(detuning)

        % Equal electron couplings produce opposite detunings on the two branches
        inter.coupling.matrix{1,3}=diag([0 0 detuning(k)]);
        inter.coupling.matrix{2,3}=diag([0 0 detuning(k)]);
        spin_system=create(sys,inter);
        echo=clexpn(spin_system,@hahn_echo,parameters,'labframe');
        loss(n,k)=1-abs(echo);

    end

end

% Reversing C rephases one proton state; reversing E interchanges the branches
loss=[flipud(loss(2:end,:));loss];
loss=[fliplr(loss(:,2:end)) loss];
axes_hz.flip=[-flip(end:-1:2) flip];
axes_hz.detuning=[-detuning(end:-1:2) detuning];

% Display the pair contribution to triplet decoherence
kfigure();
surf(axes_hz.detuning,axes_hz.flip,loss);
kxlabel('hyperfine detuning, Hz');
kylabel('flip-flop coupling, Hz');
kzlabel('pair decoherence');
kgrid();

end


