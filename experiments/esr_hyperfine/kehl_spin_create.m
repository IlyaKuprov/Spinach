% creates the Map defining the spinsystem
% is directly called for the set up
%
% input parameters:
% S: electron spin quantum number
% g: g tensor (diagonal values)
% N_SpinSys: number of individual spin systems
% N_Nuclei: number of nuclei
% Nuclei: nuclei in the spin system (e.g. '1H')
% A_ten: diagonal values of A tensor
% A_ang: Euler angles for A tensor (alpha, beta, gamma)
% Q_ten: diagonal values of Q tensor
% Q_ang: Euler angles for Q tensor (alpha, beta, gamma)
%
% output parameters:
% SpinSys: the Map containing the spinsystem parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)


function SpinSys=kehl_spin_create(S,g,N_SpinSys,N_Nuclei,Nuclei,A_ten,A_ang,Q_ten,Q_ang)
if nargin<8
    Q_ten=[];
end
if nargin<9
    Q_ang=[];
end


    % Check consistency
    grumble(S,g,N_SpinSys,N_Nuclei,Nuclei,A_ten,A_ang,Q_ten,Q_ang);
    SpinSys=containers.Map;

    SpinSys("S")=S;
    SpinSys("g")=diag(g);
    SpinSys("g_iso")=(g(1)+g(2)+g(3))/3;
    SpinSys("Ni_ENDOR")=N_Nuclei;
    SpinSys("Nuclei")=Nuclei;
    SpinSys("N_SpinSys")=N_SpinSys;

    I=zeros(N_Nuclei,1);

    for i=1:N_Nuclei
        [~,multiplicity]=spin(kehl_spin_label(Nuclei{i}));
        I(i)=(multiplicity-1)/2;
    end
    SpinSys("I")=I;

    L=zeros(N_Nuclei,12);

    for i=1:N_Nuclei
        L(i,1:3)=A_ten(3*i-2:3*i)*10^6;
        L(i,4:6)=A_ang(3*i-2:3*i);
        if size(Q_ten,2)>0
            SpinSys("Q_used")=true;
            L(i,7:9)=Q_ten(3*i-2:3*i)*10^6;
            L(i,10:12)=Q_ang(3*i-2:3*i);
        else
            SpinSys("Q_used")=false;
        end
    end


    % transform L to g-frame A and Q tensors
    [A,Q]=kehl_to_g_frame(L,N_Nuclei);
    SpinSys("A")=A;
    SpinSys("Q")=Q;
    SpinSys("EPR_Nucs_used")=false;
    SpinSys("CS_used")=false;
    SpinSys("D_used")=false;
end

function grumble(S,g,N_SpinSys,N_Nuclei,Nuclei,A_ten,A_ang,Q_ten,Q_ang)
if ~isnumeric(S)
    error('S must be numeric.');
end
if ~isnumeric(g)
    error('g must be numeric.');
end
if ~isnumeric(N_SpinSys)
    error('N_SpinSys must be numeric.');
end
if ~isnumeric(N_Nuclei)
    error('N_Nuclei must be numeric.');
end
if (~iscell(Nuclei))&&(~ischar(Nuclei))&&(~isstring(Nuclei))
    error('Nuclei must be a cell array, a character string, or a string scalar.');
end
if ~isnumeric(A_ten)
    error('A_ten must be numeric.');
end
if ~isnumeric(A_ang)
    error('A_ang must be numeric.');
end
if ~isnumeric(Q_ten)
    error('Q_ten must be numeric.');
end
if ~isnumeric(Q_ang)
    error('Q_ang must be numeric.');
end
end

