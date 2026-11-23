% Proton-pair contribution to the Hahn echo of two electron spins:
% Fig. 3(b) of Chen et al. (https://doi.org/10.1063/5.0303171).
% The two polarised triplets select the branches in Eq. (29). In the
% antiparallel proton subspace the flip-flop matrix element is C and
% the diagonal energy difference is D+E or D-E. Common Zeeman and
% exchange energies do not enter this pair echo.
% The plotted loss is 1-|L| at a total echo time of 40 microseconds.
%
% ilya.kuprov@weizmann.ac.il

function gcce_two_centres()

% Electron pair and antiparallel proton pair in the secular Hamiltonian
sys.isotopes={'E','E','1H','1H'};
sys.magnet=0.01;
sys.output='hush';
sys.tols.cce_mem=2^24;
inter.zeeman.scalar={0,0,hz2ppm(5000,sys.magnet,'1H'),0};
inter.coupling.matrix=cell(4);

% Refocus the double-quantum coherence between the two polarised triplets
parameters.observed=[1 2];
parameters.transition=struct('spins',[1 2],'ms',[0.5 -0.5;0.5 -0.5]);
parameters.times=40e-6;
parameters.cce.order=2;
parameters.cce.approximation='generalized';

% Antiparallel proton pair
parameters.bath_state={diag([1 0]),diag([0 1])};

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
flip=[-flip(end:-1:2) flip];
detuning=[-detuning(end:-1:2) detuning];

% Display the pair contribution to triplet decoherence
kfigure();
surf(detuning,flip,loss,'FaceColor',[1 0.65 0.15],'EdgeColor',[0.3 0.25 0.18]);
xlim([-1e5 1e5]);
ylim([-4000 4000]);
zlim([0 1]);
xticks(-1e5:5e4:1e5);
yticks(-4000:2000:4000);
zticks(0:0.5:1);
pbaspect([1 1 0.5]);
set(gca,'YDir','reverse','Projection','perspective');
view(-21,35);
kxlabel('$E_{-1,k}$ / Hz');
kylabel('$C_{\alpha k}$ / Hz');
kzlabel('$g_k$');
box on;
grid off;

end


