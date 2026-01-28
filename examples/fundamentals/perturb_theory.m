% Rayleigh-Schrodinger and Van Vleck perturbation theory
% modules test. Eigenvector representations differ in the
% two theories, but the energies are the same.
%
% ilya.kuprov@weizmann.ac.il

function perturb_theory()

% Settings
ham_dim=512;   % Matrix dimension
max_ord=5;     % Max RSPT order
per_amp=1/25;  % Relative H1 ampl

% H0 - Zeeman interaction
sigma=pauli(ham_dim);
H0=full(sigma.z);

% H1 - random matrix
H1=randn(ham_dim)+1i*randn(ham_dim);
H1=per_amp*(H1+H1')/2;

% Energies - perturbation theories
E_rs=zeros(ham_dim,max_ord);
E_vv=zeros(ham_dim,max_ord);
for n=1:max_ord
    E_rs(:,n)=rspert(diag(H0),H1,n);
    E_rs(:,n)=sort(E_rs(:,n),'descend');
    E_vv(:,n)=vvpert(diag(H0),H1,n);
    E_vv(:,n)=sort(E_rs(:,n),'descend');
end

% Energies - diagonalisation
E_inf=eig(H0+H1,'vector');
E_inf=sort(real(E_inf),'descend');

% Comparison with exact diagonalisation
diffs_rs=[diag(H0) real(E_rs)]-E_inf;
diffs_rs=sqrt(sum(diffs_rs.^2,1))/norm(E_inf,2);
diffs_vv=[diag(H0) real(E_vv)]-E_inf;
diffs_vv=sqrt(sum(diffs_vv.^2,1))/norm(E_inf,2);
kfigure(); subplot(1,2,1); hold on;
plot(0:max_ord,diffs_rs,'-'); 
plot(0:max_ord,diffs_vv,'--');
set(gca,'YScale','log'); kgrid;
kxlabel('perturbation theory order'); 
kylabel('energies, $\|pt-exact\|/\|exact\|$');
legend({'RSPT','VVPT'}); box on;

% Eigensystems, PTs
V_rs=cell(max_ord,1);
V_vv=cell(max_ord,1);
for n=1:max_ord

    % RSPT gets the eigensystem directly
    [~,V_rs{n}]=rspert(diag(H0),H1,n);

    % VVPT returns a generator that needs exponentiation
    [~,V_vv{n}]=vvpert(diag(H0),H1,n); V_vv{n}=expm(V_vv{n});

end

% Zero order is unit matrices
V_rs=[{eye(ham_dim)}; V_rs];
V_vv=[{eye(ham_dim)}; V_vv];

% Eigensystem - diagonalisation
[V_inf,E_inf]=eig(H0+H1,'vector');
[~,idx]=sort(real(E_inf),'descend');
V_inf=V_inf(:,idx); % Sorting

% Comparison with exact diagonalisation
diffs_rs=cellfun(@(V)norm(abs(V'*V_inf)-eye(size(V)),2),V_rs);
diffs_vv=cellfun(@(V)norm(abs(V'*V_inf)-eye(size(V)),2),V_vv);
subplot(1,2,2); hold on;
plot(0:max_ord,diffs_rs','-'); 
plot(0:max_ord,diffs_vv','--');
set(gca,'YScale','log'); kgrid;
kxlabel('perturbation theory order'); 
kylabel('eigensystem overlap residual');
legend({'RSPT','VVPT'}); box on;

end

