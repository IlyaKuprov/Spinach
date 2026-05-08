% adds nuclear dipolar coupling values to the spinsystem
% is directly called for the set up
%
% input parameters:
% spinSys_in: the Map with the spin parameters
% D_ten: the diagonal values of the coupling tensor in MHz
% D_ang: the Euler angles for the coupling tensor (alpha, beta, gamma)
%
% output parameters:
% SpinSys_out: updated Map with spin parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)


function SpinSys_out=kehl_spin_d(spinSys_in,D_ten,D_ang)

    % Check consistency
    grumble(spinSys_in,D_ten,D_ang);
    SpinSys_out=spinSys_in;
    Ni_Dip=size(D_ten,2)/3;
    SpinSys_out('D_used')=true;

    L=zeros(Ni_Dip,12);
    for i=1:Ni_Dip
            L(i,1:3)=D_ten(3*i-2:3*i)*10^3;
            L(i,4:6)=D_ang(3*i-2:3*i);
    end

    [D,~]=kehl_to_g_frame(L,Ni_Dip);

    SpinSys_out("D")=D;
end

function grumble(spinSys_in,D_ten,D_ang)
if ~isa(spinSys_in,'containers.Map')
    error('spinSys_in must be a containers.Map object.');
end
if ~isnumeric(D_ten)
    error('D_ten must be numeric.');
end
if ~isnumeric(D_ang)
    error('D_ang must be numeric.');
end
end

