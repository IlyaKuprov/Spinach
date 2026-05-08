% adds chemical shielding values to the spinsystem
% is directly called for the set up
%
% input parameters:
% spinSys_in: the Map with the spin parameters
% CS_ten: the diagonal values of the shielding tensor in ppm
% CS_ang: the Euler angles for the shielding tensor (alpha, beta, gamma)
%
% output parameters:
% SpinSys_out: updated Map with spin parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)


function SpinSys_out=kehl_spin_cs(spinSys_in,CS_ten,CS_ang)

    % Check consistency
    grumble(spinSys_in,CS_ten,CS_ang);
    SpinSys_out=spinSys_in;
    Ni_ENDOR=spinSys_in("Ni_ENDOR");
    SpinSys_out('CS_used')=true;

    L=zeros(Ni_ENDOR,12);
    for i=1:Ni_ENDOR
            L(i,1:3)=CS_ten(3*i-2:3*i)*10^6;
            L(i,4:6)=CS_ang(3*i-2:3*i);
    end

    [CS,~]=kehl_to_g_frame(L,Ni_ENDOR);

    SpinSys_out("CS")=CS;
end

function grumble(spinSys_in,CS_ten,CS_ang)
if ~isa(spinSys_in,'containers.Map')
    error('spinSys_in must be a containers.Map object.');
end
if ~isnumeric(CS_ten)
    error('CS_ten must be numeric.');
end
if ~isnumeric(CS_ang)
    error('CS_ang must be numeric.');
end
end

