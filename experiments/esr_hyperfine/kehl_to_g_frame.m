% transforms A and Q tensor into the g diagonal frame using defined Euler angles
%
% input parameters:
% L: matrix with diagonal values (i,1:3) and Euler angles (i,4:6) of nucleus i
% Ni: number of nuclei
%
% output parameters:
% A: HFC tensor in g-frame
% Q: NQC tensor in g-frame
%
% February 2024 A. Kehl (akehl@gwdg.de)


function [A,Q]=kehl_to_g_frame(L,Ni)

    % Check consistency
    grumble(L,Ni);
R=zeros(3,3);
RQ=zeros(3,3);

for k=1:Ni
    alpha=L(k,4)*pi/180;
    beta=L(k,5)*pi/180;
    gamm=L(k,6)*pi/180;
    R(1,1)=cos(beta)*cos(alpha)*cos(gamm)-sin(alpha)*sin(gamm);
    R(1,2)=cos(beta)*sin(alpha)*cos(gamm)+cos(alpha)*sin(gamm);
    R(1,3)=-sin(beta)*cos(gamm);
    R(2,1)=-cos(beta)*cos(alpha)*sin(gamm)-sin(alpha)*cos(gamm);
    R(2,2)=-cos(beta)*sin(alpha)*sin(gamm)+cos(alpha)*cos(gamm);
    R(2,3)=sin(beta)*sin(gamm);
    R(3,1)=sin(beta)*cos(alpha);
    R(3,2)=sin(beta)*sin(alpha);
    R(3,3)=cos(beta);

    X=R*diag(L(k,1:3))*R';
    A((k-1)*3+1:(k-1)*3+3,:)=X;

    if not(L(k,7:9)==[0,0,0])
        alphaQ=L(k,10)*pi/180;
        betaQ=L(k,11)*pi/180;
        gammQ=L(k,12)*pi/180;
        RQ(1,1)=cos(betaQ)*cos(alphaQ)*cos(gammQ)-sin(alphaQ)*sin(gammQ);
        RQ(1,2)=cos(betaQ)*sin(alphaQ)*cos(gammQ)+cos(alphaQ)*sin(gammQ);
        RQ(1,3)=-sin(betaQ)*cos(gammQ);
        RQ(2,1)=-cos(betaQ)*cos(alphaQ)*sin(gammQ)-sin(alphaQ)*cos(gammQ);
        RQ(2,2)=-cos(betaQ)*sin(alphaQ)*sin(gammQ)+cos(alphaQ)*cos(gammQ);
        RQ(2,3)=sin(betaQ)*sin(gammQ);
        RQ(3,1)=sin(betaQ)*cos(alphaQ);
        RQ(3,2)=sin(betaQ)*sin(alphaQ);
        RQ(3,3)=cos(betaQ);
        Y=RQ*diag(L(k,7:9))*RQ';
        Q((k-1)*3+1:(k-1)*3+3,:)=Y;
    else
        Q ((k-1)*3+1:(k-1)*3+3,:)=zeros(3,3);
    end
end
end

function grumble(L,Ni)
if ~isnumeric(L)
    error('L must be numeric.');
end
if ~isnumeric(Ni)
    error('Ni must be numeric.');
end
end

