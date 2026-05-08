% generate diagonal spin operators
%
% input parameters:
% spinOps: the Map containing the spin operators
% N_Nuclei: number of nuclei
%
% output parameters:
% Ops: the Map containing the diagonal spin operators
%
% February 2024 A. Kehl (akehl@gwdg.de)


function Ops=kehl_diag_ops(spinOps,N_Nuclei)

    % Check consistency
    grumble(spinOps,N_Nuclei);
    Ops=containers.Map;
    D=eye(size(spinOps('Sz')));


    Ops('Sx')=D\spinOps('Sx')*D;
    Ops('Sy')=D\spinOps('Sy')*D;
    Ops('Sz')=D\spinOps('Sz')*D;

    Ix=spinOps('Ix');
    Iy=spinOps('Iy');
    Iz=spinOps('Iz');

    Ix_D{N_Nuclei}=[];
    Iy_D{N_Nuclei}=[];
    Iz_D{N_Nuclei}=[];

    for i=1:size(Ix,2)
        Ix_D{i}=D/Ix{i}*D;
        Iy_D{i}=D/Iy{i}*D;
        Iz_D{i}=D/Iz{i}*D;
    end

    Ops('Ix')=Ix_D;
    Ops('Iy')=Iy_D;
    Ops('Iz')=Iz_D;
end

function grumble(spinOps,N_Nuclei)
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if ~isnumeric(N_Nuclei)
    error('N_Nuclei must be numeric.');
end
end

