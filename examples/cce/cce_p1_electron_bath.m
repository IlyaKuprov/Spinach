% P1 pair energies and resonant electron flip-flops at 500 G along [111],
% from Park et al., Fig. 5a,c,d (doi:10.1038/s41534-022-00605-4).
% The level diagrams use the diagonal high-field Zeeman and hyperfine
% terms of the product-state picture, without quadrupolar splittings
% or a common energy offset. The diagrams use the Fig. 5b nuclear state
% numbering and do not calculate the bath ensemble of Fig. 4.
%
% ilya.kuprov@weizmann.ac.il

function cce_p1_electron_bath()

% NV frame and the aligned and tilted Jahn-Teller axes
magnet=0.05;
z_axis=[1 1 1]/sqrt(3);
x_axis=[1 -1 0]/sqrt(2);
frame=[x_axis;cross(z_axis,x_axis);z_axis];
jt_axes=[1 1 1;1 -1 -1]*frame.'/sqrt(3);

% Experimental hyperfine tensors from Supplementary Note 2
hyperfine=cell(1,2);
for n=1:2

    % Rotate each axial tensor into the field frame
    axis_proj=jt_axes(n,:).'*jt_axes(n,:);
    hyperfine{n}=81.312*eye(3)+(114.0264-81.312)*axis_proj;

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
model_names={'(a) Hyperfine off','(c) Aligned JT axes','(d) One tilted JT axis'};
omega_n=-spin('14N')*magnet/(2*pi*1e6);
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

end


