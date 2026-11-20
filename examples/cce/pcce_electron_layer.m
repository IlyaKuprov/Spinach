% Partitioned Hahn echo of an NV centre and 20 bath electrons on a 20 nm
% lattice, following Schätzle et al., doi:10.1103/PhysRevB.110.L220302,
% Fig. S4(a). The finite lattice boundary, pair assignment, and random
% draws were not deposited; this reconstruction uses the nearest square
% lattice shells, one minimum-distance pairing, and equal-size constrained
% k-means partitions for the higher orders, without fitting.
% Syntax:
%
%             [answer,diagnostics]=pcce_electron_layer()
%
% Columns are pCCE(2,1), internally averaged pCCE(2,2), and normally
% averaged pCCE(2,2), pCCE(2,4), and pCCE(2,5) magnetisations. Diagnostics
% contains the echo times in seconds, lattice sites in nm, the chosen
% dimer partition, and the partitions of size four and five.
%
% Calculation time: about two minutes on a Xeon Platinum 8270.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=pcce_electron_layer()

% Place the NV at the origin with the magnetic field normal to the layer
[x,y]=ndgrid(-2:2);
keep=(x.^2+y.^2>0)&(x.^2+y.^2<=5);
coords=20*[x(keep) y(keep) zeros(20,1)];
times=(0:10:250)'*1e-6;

% Use a minimum-distance dimer partition and the authors' constrained k-means for K=4 and K=5
pairs={[1 4],[2 6],[3 8],[5 10],[7 11],[9 13],[12 17],[14 18],[15 19],[16 20]};
groups_four={[1 4 5 10],[2 3 6 7],[8 11 12 17],[9 13 14 18],[15 16 19 20]};
groups_five={[1 4 5 9 10],[2 3 6 7 8],[11 12 16 17 20],[13 14 15 18 19]};

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
partitions={num2cell(1:20),pairs,groups_four,groups_five};

% Freeze exterior bath spins along the field, but trace over retained spins
stream=RandStream('mt19937ar','Seed',1);
draws=double(rand(stream,20,100)>=0.5)-0.5;
answer=zeros(numel(times),numel(partitions)+1);
for k=1:numel(partitions)

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
    elseif k==2
        answer(:,2:3)=[internal mean(samples,2)];
    else
        answer(:,k+1)=mean(samples,2);
    end

end

% Compare partition sizes and averaging prescriptions on the total echo-time axis
diagnostics=struct('times',times,'coordinates_nm',coords,'pairs',{pairs},...
                   'groups_four',{groups_four},'groups_five',{groups_five});
kfigure();
plot(1e6*times,answer,'-o');
kxlabel('Total echo time $2\tau$ / $\mu$s');
kylabel('Normalised NV magnetisation');
ylim([0 1]);
klegend({'pCCE(2,1)','pCCE(2,2), internal average','pCCE(2,2)',...
         'pCCE(2,4)','pCCE(2,5)'});
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
H=sparse(dim,dim);
for n=1:count
    for k=1:n-1
        hopping=raising{n}*raising{k}';
        ising=spdiags(zops(:,n).*zops(:,k),0,dim,dim);
        H=H+couplings(sites(n),sites(k))*(ising-(hopping+hopping')/4);
    end
end

% Only exterior spins generate frozen fields; NV projections are 0 and -1
outside=setdiff(1:20,sites);
offsets=zops*couplings(sites,outside)*draws(outside,:);
shift=zops*hyperfine(sites);

% The two fully polarised sectors each contribute one to the trace
traces(:)=2;
magnetisation=sum(zops,2);
for m=1-count/2:count/2-1

    % Both branch Hamiltonians conserve bath magnetisation
    sector=(magnetisation==m);
    block=full(H(sector,sector));
    nv_shift=diag(shift(sector));
    for n=1:size(draws,2)

        % Diagonalise both branches and express the second in the first eigenbasis
        branch=block+diag(offsets(sector,n));
        [va,ea]=eig(branch,'vector');
        [vb,eb]=eig(branch-nv_shift,'vector');
        overlap=va'*vb;
        gaps=ea-ea.';
        phase=(times*eb.')/2;

        % When the first propagator is diagonal, only second-propagator intensities enter
        for k=2:numel(times)
            real_part=(overlap.*cos(phase(k,:)))*overlap';
            imag_part=(overlap.*sin(phase(k,:)))*overlap';
            intensities=real_part.^2+imag_part.^2;
            traces(k,n)=traces(k,n)+sum(cos((times(k)/2)*gaps).*intensities,'all');
        end

    end

end

% Normalise by the full Hilbert-space dimension, including the time-zero identity
traces=traces/dim;
traces(1,:)=1;

end


