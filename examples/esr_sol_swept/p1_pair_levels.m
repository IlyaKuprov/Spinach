% P1 pair energies and resonant electron flip-flops at 500 G along [111],
% from Park et al., Fig. 5a,c,d (doi:10.1038/s41534-022-00605-4).
% The level diagrams use the diagonal high-field Zeeman and hyperfine
% terms of the product-state picture, without quadrupolar splittings
% or a common energy offset. The diagrams use the Fig. 5b nuclear state
% numbering and do not calculate the bath ensemble of Fig. 4.
%
% ilya.kuprov@weizmann.ac.il

function p1_pair_levels()

% NV frame and the aligned and tilted Jahn-Teller axes
sys.magnet=0.05;
sys.isotopes={'E','E','14N','14N'};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
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

% Fig. 5b nuclear ordering within the two opposite electron configurations
projections=[-1 -1; 0 -1; -1 0; 1 -1; 0 0; -1 1; 1 0; 0 1; 1 1];
nuc_idx=3*(1-projections(:,1))+2-projections(:,2);
states=[9+nuc_idx 18+nuc_idx];

% Model choices in panel order a, c, d
axis_ids=[1 1;1 1;1 2];
hf_on=[false true true];
model_names={'(a) Hyperfine off','(c) Aligned JT axes','(d) One tilted JT axis'};
energies=zeros(9,2,3);
channels=false(9,9,3);
for n=1:3

    % Set each electron-nitrogen hyperfine tensor in Hz
    inter.coupling.matrix=cell(4);
    for k=1:2
        inter.coupling.matrix{k,k+2}=hf_on(n)*1e6*hyperfine{axis_ids(n,k)};
    end
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    spin_system=assume(spin_system,'labframe');

    % Read diabatic product-state energies, not mixed eigenvalues
    [H,Q]=hamiltonian(spin_system);
    diagonal=full(diag(H+orientation(Q,[0 0 0])))/(2*pi*1e6);
    levels=diagonal(states);

    % Electron flip-flops conserve both nuclear projections
    F=operator(spin_system,{'L+','L-'},{1,2});
    allowed=full(F(states(:,1),states(:,2))~=0);
    energies(:,:,n)=levels;

    % Cancellation of the common electron Zeeman terms sets the roundoff scale
    tolerance=size(H,1)*eps(norm(diagonal,inf));
    channels(:,:,n)=allowed&(abs(levels(:,1)-levels(:,2).')<=tolerance);

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


