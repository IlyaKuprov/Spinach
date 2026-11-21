% Projected CCE-2 Hahn echo for the proton-only middle panel of
% Ryan et al., Figure 3, doi:10.1063/5.0278659. The selected levels
% are the 4th and 12th ascending eigenstates of the electron–
% vanadium-51 register at 0.33 T.
% Register-changing bath matrix elements are omitted, but all
% conditional hyperfine and proton-proton dipolar terms are kept.
% The bath is unpolarised, and the ideal pulse swaps the two levels.
% The linked cluster factors use the full complex Hahn overlap:
% L=prod(L_j)*prod(L_jk/(L_j*L_k)) for singles j and pairs jk.
%
%                     [echo,times]=votpp_proton_echo()
%
% Row vectors: normalised complex coherence and total echo time
% in seconds. Internally the Hamiltonians are in radians per ms.
% Crystal sites and tensors are from the LunghiGroup companion
% archive Spin-Decoherence---CCE (VOTPP); no trace data are used.

function [echo,times]=votpp_proton_echo()

% Physical parameters in gauss, angstroms, and radian-millisecond units
field=3300;
gyro_h=26.7522;
gyro_v=-7.05;
gyro_e=-8794.10005384623*[ 1.9677654  0.0000435 -0.0000131;
                          0.0000434  1.9838165  0.0000203;
                         -0.0000436 -0.0000195  1.9838159];
A=2*pi*1e3*[-472.6908    0.0002    -0.0001;
               0.0002 -165.5146     0.0005;
              -0.0001   -0.0005  -165.5145];
quad=-350*2*pi;

% Assemble the coupled register and project its spin operators
proton=pauli(2);
vanadium=pauli(8);
proton_ops=cat(3,full(proton.x),full(proton.y),full(proton.z));
labels='xyz';
nuc_ops=zeros(16,16,3);
ele_ops=zeros(16,16,3);
for k=1:3
    nuc_ops(:,:,k)=kron(full(vanadium.(labels(k))),eye(2));
    ele_ops(:,:,k)=kron(eye(8),proton_ops(:,:,k));
end
H=-gyro_v*field*nuc_ops(:,:,3)+quad*(nuc_ops(:,:,3)^2-3.5*4.5*eye(16)/3);
for k=1:3

    % Add electronic Zeeman and electron-vanadium hyperfine couplings
    H=H+field*gyro_e(3,k)*ele_ops(:,:,k);
    for n=1:3
        H=H+A(n,k)*nuc_ops(:,:,n)*ele_ops(:,:,k);
    end

end
[states,energies]=eig(H,'vector');
[~,order]=sort(energies,'ascend');
states=states(:,order([4 12]));
moments=zeros(2,3);
for k=1:3

    % Gyromagnetically weighted spin expectations of the two register levels
    for n=1:3
        moments(:,k)=moments(:,k)+gyro_e(n,k)*real(diag(states'*ele_ops(:,:,n)*states));
    end
    moments(:,k)=moments(:,k)+gyro_v*real(diag(states'*nuc_ops(:,:,k)*states));

end

% Translate the 56 relative proton sites and select the 20 angstrom bath
sites=readmatrix(fullfile(fileparts(mfilename('fullpath')),'votpp_proton_sites.csv'));
cell_vec=[9.67194       0.0          0.0;
          0.000288633  13.2521       0.0;
          0.00231045   -0.0000010487 13.2521];
bath_radius=20;
cell_inv=inv(cell_vec);
limits=ceil(bath_radius*vecnorm(cell_inv,2,2).'+max(abs(sites*cell_inv.'),[],1));
[x,y,z]=ndgrid(-limits(1):limits(1),-limits(2):limits(2),-limits(3):limits(3));
translations=[x(:) y(:) z(:)]*cell_vec.';
coords=reshape(reshape(sites,[],1,3)+reshape(translations,1,[],3),[],3);
coords=coords(vecnorm(coords,2,2)<bath_radius,:);

% Select the proton pairs within 8 angstroms
separation=reshape(coords,[],1,3)-reshape(coords,1,[],3);
[left,right]=find(triu(sum(separation.^2,3)<8^2,1));

% Construct the full conditional one- and two-proton Hamiltonians
hf=pagemtimes(moments,dipole_tensor(coords))*gyro_h;
hf(:,3,:)=hf(:,3,:)-gyro_h*field;
dip=dipole_tensor(coords(left,:)-coords(right,:))*gyro_h^2;
pair_coupling=zeros(4,4,numel(left));
for k=1:3
    for n=1:3
        pair_coupling=pair_coupling+kron(proton_ops(:,:,k),proton_ops(:,:,n)).*dip(k,n,:);
    end
end
singles=cell(1,2);
pairs=cell(1,2);
for k=1:2

    % Embed the two local Hamiltonians and the shared dipolar coupling
    singles{k}=sum(proton_ops.*reshape(hf(k,:,:),1,1,3,[]),3);
    singles{k}=reshape(singles{k},2,2,[]);
    pairs{k}=pair_coupling;
    for n=1:3
        pairs{k}=pairs{k}+kron(proton_ops(:,:,n),eye(2)).*hf(k,n,left)+...
                          kron(eye(2),proton_ops(:,:,n)).*hf(k,n,right);
    end

end

% Hahn overlaps and the linked CCE-2 product on the figure time axis
times=linspace(0,30e-6,121);
single_echo=cluster_echo(singles{1},singles{2},times*1e3);
pair_echo=cluster_echo(pairs{1},pairs{2},times*1e3);
echo=prod(single_echo,1).*prod(pair_echo./(single_echo(left,:).*single_echo(right,:)),1);

end

% Magnetic point-dipole tensor before gyromagnetic factors
function T=dipole_tensor(coords)

% The prefactor is hbar times mu0/(4*pi) in the internal units
radius=vecnorm(coords,2,2);
direction=reshape((coords./radius).',3,1,[]);
T=1.05457172*(eye(3)-3*pagemtimes(direction,'none',direction,'transpose'))./...
  reshape(radius.^3,1,1,[]);

end

% Compute normalised Hahn overlaps for pages of conditional Hamiltonians
function echo=cluster_echo(Ha,Hb,times)

% Diagonalise the two conditional Hamiltonians once per cluster
[Va,Ea]=pageeig(Ha,'vector');
[Vb,Eb]=pageeig(Hb,'vector');
echo=zeros(size(Ha,3),numel(times));
for k=1:numel(times)

    % Two bath histories spend the same half-interval on opposite levels
    Ua=pagemtimes(Va.*pagetranspose(exp(-0.5i*times(k)*Ea)),'none',Va,'ctranspose');
    Ub=pagemtimes(Vb.*pagetranspose(exp(-0.5i*times(k)*Eb)),'none',Vb,'ctranspose');
    ab=pagemtimes(Ua,Ub);
    ba=pagemtimes(Ub,Ua);

    % Use the squared norm of the unitary difference to avoid cancellation
    overlap=1-sum(abs(ab-ba).^2,[1 2])/(2*size(Ha,1))+...
            1i*sum(imag(ab.*conj(ba)),[1 2])/size(Ha,1);
    echo(:,k)=overlap(:);

end

end


