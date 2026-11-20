% Partitioned Hahn echo of an NV centre and 20 bath electrons on a 20 nm
% lattice, following Schätzle et al., doi:10.1103/PhysRevB.110.L220302,
% Fig. S4(a). The finite lattice boundary, pair assignment, and random
% draws were not deposited; this reconstruction uses the nearest square
% lattice shells and one minimum-distance pairing, without fitting.
% Syntax:
%
%             [answer,diagnostics]=pcce_electron_layer()
%
% Columns are pCCE(2,1), internally averaged pCCE(2,2), and normally
% averaged pCCE(2,2) magnetisations. Diagnostics contains the echo times
% in seconds, the lattice sites in nm, and the chosen dimer partition.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=pcce_electron_layer()

% Place the NV at the origin with the magnetic field normal to the layer
[x,y]=ndgrid(-2:2);
keep=(x.^2+y.^2>0)&(x.^2+y.^2<=5);
coords=20*[x(keep) y(keep) zeros(20,1)];
pairs={[1 4],[2 6],[3 8],[5 10],[7 11],[9 13],[12 17],[14 18],[15 19],[16 20]};
times=(0:10:250)'*1e-6;

% Secular dipolar constants in rad/s for the perpendicular field
couplings=zeros(20);
hyperfine=zeros(20,1);
for n=1:20
    hyperfine(n)=xyz2dd([0 0 0],10*coords(n,:),'E','E');
    for k=1:n-1
        couplings(k,n)=xyz2dd(10*coords(k,:),10*coords(n,:),'E','E');
        couplings(n,k)=couplings(k,n);
    end
end

% Enumerate a complete group graph with the NV distinct from the bath
model=struct('physical_ids',1:21,'observed_ids',21,'bath_ids',1:20,...
             'multiplicities',[2*ones(1,20) 3]);
graph=struct('bath_ids',1:20,'adjacency',~eye(20),'scores',couplings);
selection=struct('order',2,'expansion','partition');
partitions={num2cell(1:20),pairs};

% Freeze exterior bath spins along the field, but trace over retained spins
stream=RandStream('mt19937ar','Seed',1);
draws=double(rand(stream,20,100)>=0.5)-0.5;
answer=zeros(numel(times),3);
for k=1:2

    % Generate all one- and two-partition clusters and their Möbius powers
    selection.groups=partitions{k};
    clusters=cce_clusters(model,graph,selection);
    samples=ones(numel(times),size(draws,2));
    internal=ones(numel(times),1);

    % Average raw traces before reduction only for the internal-average curve
    for n=1:numel(clusters)
        sites=clusters(n).sites;
        traces=cluster_echo(sites,couplings,hyperfine,draws,times);
        power=double(clusters(n).mobius_power);
        samples=samples.*traces.^power;
        if k==2
            internal=internal.*mean(traces,2).^power;
        end
    end
    if k==1
        answer(:,1)=mean(samples,2);
    else
        answer(:,2:3)=[internal mean(samples,2)];
    end

end

% Compare the two averaging prescriptions on the total echo-time axis
diagnostics=struct('times',times,'coordinates_nm',coords,'pairs',{pairs});
kfigure();
plot(1e6*times,answer,'-o');
kxlabel('Total echo time $2\tau$ / $\mu$s');
kylabel('Normalised NV magnetisation');
ylim([0 1]);
klegend({'pCCE(2,1)','pCCE(2,2), internal average','pCCE(2,2)'});
kgrid();

end

% Infinite-temperature trace of the two conditional Hahn-echo paths
function traces=cluster_echo(sites,couplings,hyperfine,draws,times)

% Singleton longitudinal Hamiltonians commute on both NV branches
count=numel(sites);
traces=ones(numel(times),size(draws,2));
if count==1
    return
end

% Embed the retained spin operators in the local Hilbert space
ops=pauli(2);
dim=2^count;
zops=zeros(dim,count);
raising=cell(1,count);
for n=1:count
    left=speye(2^(n-1));
    right=speye(2^(count-n));
    zops(:,n)=diag(kron(kron(left,ops.z),right));
    raising{n}=kron(kron(left,ops.p),right);
end

% Retain the secular Ising and homonuclear flip-flop terms of Eq. S7
H=zeros(dim);
for n=1:count
    for k=1:n-1
        hopping=raising{n}*raising{k}';
        H=H+couplings(sites(n),sites(k))*(diag(zops(:,n).*zops(:,k))-(hopping+hopping')/4);
    end
end

% Only exterior spins generate frozen fields; NV projections are 0 and -1
outside=setdiff(1:20,sites);
offsets=zops*couplings(sites,outside)*draws(outside,:);
shift=diag(zops*hyperfine(sites));
for n=1:size(draws,2)

    % Diagonalise both branches once for the complete time axis
    branch=H+diag(offsets(:,n));
    [va,ea]=eig(branch,'vector');
    [vb,eb]=eig(branch-shift,'vector');
    for k=1:numel(times)
        ua=(va.*exp(-0.5i*times(k)*ea.'))*va';
        ub=(vb.*exp(-0.5i*times(k)*eb.'))*vb';
        traces(k,n)=real(sum(conj(ua*ub).*(ub*ua),'all'))/dim;
    end

end

end


