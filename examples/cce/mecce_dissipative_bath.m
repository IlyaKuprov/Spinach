% Central-spin free-induction coherence in a dissipative 12-spin XXZ ring,
% reproducing the five ME-CCE curves of Onizhuk et al., Phys. Rev. Lett.
% 132, 250401 (2024), Fig. 1(c), https://doi.org/10.1103/PhysRevLett.132.250401.
% Table S1 supplies the couplings. Fig. 1(c) is matched on its displayed time
% grid by using the tabulated frequencies directly and a rate of 0.01 per
% unit time for each of the two jumps on each site. The blue markers are
% digitised from the paper's full-bath calculation, not simulated here.
%
% ilya.kuprov@weizmann.ac.il

function mecce_dissipative_bath()

% Central-spin and periodic nearest-neighbour bath couplings from Table S1
probe_coupling=[1.023643 1.900927 0.288319 1.897299 0.623663 0.846653 ...
                1.655405 0.818398 1.099187 0.055118 1.507026 1.076287];
bath_coupling=[0.038499 0.049414 0.029822 0.073932 0.076897 0.097981 ...
               0.071064 0.050969 0.085933 0.042525 0.092109 0.074284];
nbath=numel(bath_coupling);
sites=2:nbath+1;
neighbours=[sites(2:end) sites(1)];

% Zeeman-free spin halves with longitudinal probe coupling and bath XXZ exchange
sys.magnet=0;
sys.isotopes=[{'E'} repmat({'13C'},1,nbath)];
sys.output='hush';
inter.zeeman.matrix=repmat({zeros(3)},1,nbath+1);
inter.coupling.matrix=cell(nbath+1);
for n=1:nbath

    % Isolate the two central-spin branches while allowing nearest-neighbour flips
    inter.coupling.matrix{1,sites(n)}=diag([0 0 probe_coupling(n)]);
    inter.coupling.matrix{sites(n),neighbours(n)}=bath_coupling(n)*diag([1 1 -2]);

end

% Néel product state and equal raising and lowering jump rates
bath_rho=repmat({diag([1 0]),diag([0 1])},1,nbath/2);
jumps=repmat(struct('support',0,'operator','L+','rate',0.01),1,2*nbath);
for n=1:nbath

    % Set the two local jumps on the same physical site
    jumps(2*n-1).support=sites(n);
    jumps(2*n).support=sites(n);
    jumps(2*n).operator='L-';

end

% Specify the periodic cluster graph and frozen exterior Néel state
adjacency=sparse(1:nbath,[2:nbath 1],true,nbath,nbath);
parameters.cce.graph=struct('mode','explicit','adjacency',adjacency|adjacency.');
parameters.cce.expansion='periodic_segments';
parameters.cce.order=5;
parameters.cce.bath_spins=sites;
parameters.cce.bath_state=struct('mode','explicit_product','rho',{bath_rho});
parameters.cce.meanfield='static';
parameters.cce.approximation='projected';

% Simulate the projected free induction with native retained-cluster jumps
times=0:0.1:20;
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[0.5 -0.5]);
parameters.times=times;
parameters.output='signal';
parameters.bath_jumps=jumps;

% Allow for accumulated error in the order-five product and its denominator
parameters.cce.numerics=struct('atol',1e-7,'rtol',1e-8,'max_memory_bytes',2^29);
spin_system=create(sys,inter);
[~,details]=clexpn(spin_system,@free_evolution,parameters,'labframe');
coherence=abs(reshape(details.reduction.by_order(:,1,2:end),numel(times),[]));

% Plot the five ME-CCE curves and the published exact marker centres
exact=readmatrix(fullfile(fileparts(mfilename('fullpath')),...
                          'mecce_dissipative_bath_exact.csv'));
kfigure();
lines=plot(times,coherence);
hold on;
exact_markers=plot(exact(:,1),exact(:,2),'.','Color',[0 0.65 0.78]);
kxlabel('Time ($2\pi/\langle a_j\rangle$)');
kylabel('Coherence');
xlim([0 20]);
ylim([0 1.03]);
klegend([lines;exact_markers],{'ME-CCE1','ME-CCE2','ME-CCE3',...
                             'ME-CCE4','ME-CCE5','Exact (digitised)'},...
        'Location','southwest');
kgrid();

% Compare each order to the exact marker centres at matching time coordinates
[~,exact_idx]=ismembertol(exact(:,1),times,eps(max(times)),'DataScale',1);
inset=axes('Position',[0.58 0.59 0.29 0.29]);
plot(inset,exact(:,1),coherence(exact_idx,:)-exact(:,2));
xlim(inset,[0 20]);
ylim(inset,[-0.1 0.1]);
kylabel(inset,'$\Delta$');

end


