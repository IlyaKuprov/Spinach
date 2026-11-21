% P1 pair energies and resonant electron flip-flops at 500 G along [111],
% from Park et al., Fig. 5a,c,d (doi:10.1038/s41534-022-00605-4).
% The level diagrams use the diagonal high-field Zeeman and hyperfine
% terms of the product-state picture, without quadrupolar splittings
% or a common energy offset. The separate projected NV Hahn echo uses
% the full on-site tensors and secular inter-centre dipoles for one
% illustrative pair, not the ensemble echo of Fig. 4.
%
% Syntax: [answer,diagnostics]=cce_p1_electron_bath()
%
% Columns of answer are the complex Hahn coherences for hyperfine off,
% aligned axes, and one tilted axis, respectively. Diagnostics contain
% total echo times (s), valid time masks, cluster counts, the nuclear
% projections in Fig. 5b order, and the two sets of nine energies (MHz)
% and the 9-by-9 resonance matrix for each model.
% Calculation time: about 30 seconds including local pool startup.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_p1_electron_bath()

% NV frame and the aligned and tilted Jahn-Teller axes
sys.magnet=0.05;
z_axis=[1 1 1]/sqrt(3);
x_axis=[1 -1 0]/sqrt(2);
frame=[x_axis;cross(z_axis,x_axis);z_axis];
jt_axes=[1 1 1;1 -1 -1]*frame.'/sqrt(3);

% Experimental hyperfine and DFT quadrupole tensors from Supplementary Note 2
hyperfine=cell(1,2);
quadrupole=cell(1,2);
for n=1:2

    % Rotate each axial tensor into the field frame
    axis_proj=jt_axes(n,:).'*jt_axes(n,:);
    hyperfine{n}=81.312*eye(3)+(114.0264-81.312)*axis_proj;
    quadrupole{n}=1.255*(eye(3)-3*axis_proj);

end

% Spin operators and the nuclear product basis in Fig. 5b order
S=pauli(2);
I=pauli(3);
projections=[-1 -1; 0 -1; -1 0; 1 -1; 0 0; -1 1; 1 0; 0 1; 1 1];
nuc_idx=3*(1-projections(:,1))+2-projections(:,2);
Z=full(diag(I.z));
z_nuc=[repelem(Z,3) repmat(Z,3,1)];
z_nuc=z_nuc(nuc_idx,:);
electron_idx=[2 3];
z_elec=[kron(diag(S.z),ones(2,1)) kron(ones(2,1),diag(S.z))];
z_elec=full(z_elec(electron_idx,:));
F=kron(kron(S.p,S.m)+kron(S.m,S.p),kron(I.u,I.u));
allowed=full(F(9*(electron_idx(1)-1)+nuc_idx,9*(electron_idx(2)-1)+nuc_idx)~=0);

% Model choices in panel order a, c, d
axis_ids=[1 1;1 1;1 2];
hf_on=[false true true];
model_names={'Hyperfine off','Aligned JT axes','One tilted JT axis'};
omega_n=-spin('14N')*sys.magnet/(2*pi*1e6);
energies=zeros(9,2,3);
channels=false(9,9,3);
for n=1:3

    % Diabatic pair energies and resonances to floating-point precision
    azz=hf_on(n)*cellfun(@(tensor) tensor(3,3),hyperfine(axis_ids(n,:)));
    levels=omega_n*sum(z_nuc,2)+(z_nuc.*azz)*z_elec.';
    energies(:,:,n)=levels;
    channels(:,:,n)=allowed&(abs(levels(:,1)-levels(:,2).')<=eps(max(abs(levels),[],'all')));

end

% Plot both electron configurations at each nuclear product state
kfigure('Position',[100 100 1200 380]);
tiledlayout(1,3);
for n=1:3

    % Blue segments mark resonant pairs linked by the flip-flop operator
    nexttile();
    hold on;
    [left,right]=find(channels(:,:,n));
    plot([left-0.26 right+0.26].', ...
         [energies(left,1,n) energies(right,2,n)].','b-','LineWidth',1);
    plot((1:9)-0.26,energies(:,1,n),'k_','MarkerSize',6,'LineWidth',1.5);
    plot((1:9)+0.26,energies(:,2,n),'r_','MarkerSize',6,'LineWidth',1.5);
    hold off;
    xlim([0.5 9.5]);
    ylim padded;
    xticks(1:9);
    kxlabel('Fig. 5b state index (black: unprimed; red: primed)');
    kylabel('Relative energy / MHz');
    ktitle(sprintf('%s (%d active)',model_names{n},nnz(channels(:,:,n))));

end

% One explicit diamond-site pair with co-located P1 electrons and nuclei
sys.isotopes={'E3','E','14N','E','14N'};
sys.output='hush';
sys.parallel={'processes',1};
coords=3.567*[0 0 0;6 1 1;6 1 1;5 12 2;5 12 2]*frame.';
g_e=2.00231930436256*eye(3);
inter.zeeman.matrix={g_e,g_e,zeros(3),g_e,zeros(3)};
inter.coupling.matrix=cell(5);
inter.coupling.matrix{1,1}=(2.87e9/3)*diag([-1 -1 2]);

% Retain the observed-Sz sector and high-field secular bath dipoles
for a=1:5

    % Skip co-located electron and nitrogen spins of the same centre
    for b=a+1:5

        % Convert the dipolar tensor to Hz after the secular projection
        if isequal(coords(a,:),coords(b,:))
            continue
        end
        [~,~,~,~,D]=xyz2dd(coords(a,:),coords(b,:),sys.isotopes{a},sys.isotopes{b});
        if a==1
            D(:,1:2)=0;
        elseif strcmp(sys.isotopes{a},sys.isotopes{b})
            D=D(3,3)*diag([-1/2 -1/2 1]);
        else
            D=D(3,3)*diag([0 0 1]);
        end
        inter.coupling.matrix{a,b}=D/(2*pi);

    end

end

% Project the NV transition and expand to the complete two-centre cluster
parameters.observed=1;
parameters.transition.spins=1;
parameters.transition.ms=[0 1];
parameters.times=linspace(0,100e-6,201);
parameters.cce.bath_spins=2:5;
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.expansion='partition';
parameters.cce.groups={[2 3],[4 5]};
parameters.cce.partition=struct('mode','user','max_dimension',6);
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';

% Keep the full internal P1 Hamiltonians in the echo calculation
answer=complex(zeros(numel(parameters.times),3));
valid=false(size(answer));
cluster_counts=zeros(1,3);
for n=1:3

    % Set the on-site tensors for the same pair in each model
    for k=1:2

        % Use the same quadrupole tensor with and without hyperfine coupling
        inter.coupling.matrix{2*k+1,2*k+1}=1e6*quadrupole{axis_ids(n,k)};
        inter.coupling.matrix{2*k,2*k+1}=1e6*hf_on(n)*hyperfine{axis_ids(n,k)};

    end

    % Run the native projected Hahn sequence
    spin_system=create(sys,inter);
    [answer(:,n),audit]=clexpn(spin_system,@hahn_echo,parameters,'labframe');
    valid(:,n)=audit.reduction.valid;
    cluster_counts(n)=numel(audit.clusters);

end

% Report the pair levels and the separate Hahn result
diagnostics=struct('times',parameters.times(:),'valid',valid,...
                   'cluster_counts',cluster_counts,'projections',projections,...
                   'energies',energies,'channels',channels);
kfigure();
plot(1e6*parameters.times,real(answer));
kxlabel('Total echo time / microseconds');
kylabel('Normalised Hahn-echo coherence');
ktitle('One explicit P1 pair (not Fig. 4)');
klegend(model_names{:});
kgrid();

end


