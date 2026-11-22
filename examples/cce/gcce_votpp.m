% Field-dependent energy levels of the electron–vanadium-51 register
% in VO(TPP), using Eq. (10) and Table I of doi:10.1063/5.0278659
% to reconstruct Figure 2. The field is transverse to the unique
% axis of the g, hyperfine, and quadrupolar tensors.
%
% The proton Hahn echo (middle of Figure 3) uses conditional CCE-2
% with the source crystal coordinates and tensors, retaining every
% proton within 20 Å and pairs within 8 Å. The companion input swaps
% x and z in g and A only; the vanadium quadrupole is axial along z.
% The ideal pulse swaps the 4th and 12th register levels. The bath
% is unpolarised and register-changing bath couplings are neglected.
%
% ilya.kuprov@weizmann.ac.il

function gcce_votpp()

% The unique axis is x in the crystal frame; the applied field is along z
sys.isotopes={'E','51V'};
sys.output='hush';
inter.zeeman.matrix={diag([1.968 1.984 1.984]),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([-473 -166 -166])*1e6;
inter.coupling.matrix{2,2}=-0.35e6*diag([2 -1 -1])/3;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Construct the stock laboratory-frame register Hamiltonian at the field endpoints
field_ends=[0 0.15];
H=cell(1,2);
for k=1:2

    % Build only the electron–vanadium basis
    sys.magnet=field_ends(k);
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    spin_system=assume(spin_system,'labframe');
    [iso,aniso]=hamiltonian(spin_system);
    H{k}=iso+orientation(aniso,[0 0 0]);

end

% Interpolate the Hamiltonian rather than the level energies
zeeman=(H{2}-H{1})/field_ends(2);
fields=linspace(field_ends(1),field_ends(2),100);
energies=zeros(numel(fields),size(H{1},1));
for k=1:numel(fields)

    % Sort the energy levels in spectroscopic units
    energies(k,:)=sort(eig(full(H{1}+fields(k)*zeeman))).'/(2*pi*1e9);

end

% Plot the field dependence
kfigure();
plot(fields,energies);
kxlabel('magnetic field / T');
kylabel('energy / GHz');
kgrid();

% Translate the crystal sites and retain the 20 Å proton bath
sites=readmatrix(fullfile(fileparts(mfilename('fullpath')),'votpp_proton_sites.csv'));
cell_vec=[9.67194       0.0          0.0;...
          0.000288633  13.2521       0.0;...
          0.00231045   -0.0000010487 13.2521];
bath_radius=20;
cell_inv=inv(cell_vec);
limits=ceil(bath_radius*vecnorm(cell_inv,2,2).'+max(abs(sites*cell_inv.'),[],1));
[x,y,z]=ndgrid(-limits(1):limits(1),-limits(2):limits(2),-limits(3):limits(3));
translations=[x(:) y(:) z(:)]*cell_vec.';
coords=reshape(reshape(sites,[],1,3)+reshape(translations,1,[],3),[],3);
coords=coords(vecnorm(coords,2,2)<bath_radius,:);

% Use the companion model's signed gyromagnetic tensors in rad/(s T)
gyro_h=26.7522e7;
gyro_v=-7.05e7;
gyro_e=-8794.10005384623e7*[1.9677654  0.0000435 -0.0000131;...
                           0.0000434  1.9838165  0.0000203;...
                          -0.0000436 -0.0000195  1.9838159];
sys.magnet=0.33;
sys.isotopes=[{'E','51V'} repmat({'1H'},1,size(coords,1))];

% Supply the archive Zeeman coefficients in Spinach's g and ppm units
inter.zeeman.matrix=[{gyro_e/(-spin('E')/spin_system.tols.freeg),...
                      1e6*(-gyro_v/(-spin('51V'))-1)*eye(3)} ...
                      repmat({1e6*(gyro_h/spin('1H')-1)*eye(3)},1,size(coords,1))];

% Supply the register interactions in Hz in the companion frame
nspins=numel(sys.isotopes);
inter.coupling.matrix=cell(nspins,nspins);
inter.coupling.matrix{1,2}=1e6*[-472.6908    0.0002   -0.0001;...
                                  0.0002 -165.5146    0.0005;...
                                 -0.0001   -0.0005 -165.5145];
inter.coupling.matrix{2,2}=-350e3*diag([-1 -1 2])/3;

% Build the physical dipolar tensors with both register spins at the crystal origin
pair_radius=8;
separation=reshape(coords,[],1,3)-reshape(coords,1,[],3);
[left,right]=find(triu(sum(separation.^2,3)<pair_radius^2,1));
vectors=[coords; coords(left,:)-coords(right,:)];
radius=vecnorm(vectors,2,2);
direction=reshape((vectors./radius).',3,1,[]);
tensor=(1.05457172e-11*gyro_h/(2*pi))*...
       (eye(3)-3*pagemtimes(direction,'none',direction,'transpose'))./reshape(radius.^3,1,1,[]);
for k=1:size(coords,1)

    % Coupling cells follow the stock transposed Cartesian convention
    inter.coupling.matrix{1,k+2}=(gyro_e*tensor(:,:,k)).';
    inter.coupling.matrix{2,k+2}=gyro_v*tensor(:,:,k);

end
for k=1:numel(left)

    % Supply each pair tensor once, in Hz
    inter.coupling.matrix{left(k)+2,right(k)+2}=gyro_h*tensor(:,:,size(coords,1)+k);

end

% Keep the supplied dipolar tensors rather than adding coordinate-generated copies
sys.tols.prox_cutoff=0;
inter.coordinates=num2cell([zeros(2,3); coords],2).';
spin_system=create(sys,inter);

% Evaluate the normalised complex Hahn echo with diagonal register projection
parameters.observed=[1 2];
parameters.transition=struct('spins',[1 2],'levels',[4 12]);
parameters.times=linspace(0,30e-6,121);
parameters.cce.approximation='projected';
parameters.cce.projection='conditional';
parameters.cce.order=2;
parameters.cce.parallel=struct('mode','clusters','batch_size',1+size(coords,1)+numel(left));
parameters.cce.numerics=struct('atol',1e-2,'rtol',0,'max_memory_bytes',64*2^30);
parameters.cce.graph=struct('mode','distance','radius',pair_radius);
echo=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the proton-limited coherence
kfigure();
plot(1e6*parameters.times,real(echo));
kxlabel('total echo time / $\mu$s');
kylabel('normalised coherence');
kgrid();

end


