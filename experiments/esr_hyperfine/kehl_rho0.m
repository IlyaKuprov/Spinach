% sets the initial densitiy matrix either temperature dependent or as
% approximation (-Sz)
%
% input parameters:
% constants: the Map containing the constants
% paramsENDOR: the Map containing the ENDOR parameters
% B: the main field
% geff: the effective g value
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% opt: the Map containing the optional paramters
% HF_zz, HF_zy, HF_zx: the effective HF coupling value
% NQI_zz: the effective NQ coupling value
%
% output parameters:
% rho0: the starting density matrix
%
% February 2024 A. Kehl (akehl@gwdg.de)


function [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,opt,HF_zz,HF_zy,HF_zx,NQI_zz)
if nargin<11
    NQI_zz=0;
end


    % Check consistency
    grumble(constants,paramsENDOR,B,geff,spinOps,spinSys,opt,HF_zz,HF_zy,HF_zx,NQI_zz);
 % getting parameters from Maps
 I=spinSys("I");
 Sz=spinOps('Sz');

 Ix=spinOps('Ix');
 Iy=spinOps('Iy');
 Iz=spinOps('Iz');

 v_L=paramsENDOR("v_L");
 if NQI_zz==0
     NQI_zz=zeros(spinSys("Ni_ENDOR"));
 end

 % setting the spin Hamiltonian
 H_EZ=B*geff*constants("MU_B")/constants("H")*Sz;
 H_NZ=zeros(size(H_EZ));
 H_HF=zeros(size(H_EZ));
 H_NQI=zeros(size(H_EZ));

 if spinSys('N_SpinSys')>1
    m=1;
    H_NZ=H_NZ+v_L(m)*Iz{m};
    H_HF=H_HF+HF_zz(m)*Sz*Iz{m}+HF_zy(m)*(Sz*Iy{m})+HF_zx(m)*(Sz*Ix{m});
    H_NQI=H_NQI+1/2*NQI_zz(m)*(3*Iz{m}*Iz{m}-I(m)*(I(m)+1)*eye(size(Iz{m})));
 else
    for m=1:spinSys('Ni_ENDOR')
        H_NZ=H_NZ+v_L(m)*Iz{m};
        H_HF=H_HF+HF_zz(m)*Sz*Iz{m}+HF_zy(m)*(Sz*Iy{m})+HF_zx(m)*(Sz*Ix{m});
        H_NQI=H_NQI+1/2*NQI_zz(m)*(3*Iz{m}*Iz{m}-I(m)*(I(m)+1)*eye(size(Iz{m})));
    end
 end


H_S=H_EZ+H_NZ+H_HF+H_NQI;

% calc. density matrix
if opt("temp_eff")==true

    % calculate Boltzmann factor
    Boltz=expm(-H_S/(constants("K_B")*opt("T")));

    % calculate density matrix with Temp Effect
    rho0=Boltz/trace(Boltz);
else

    % calculate density matrix without Temp Effect
    rho0=-Sz;
end

rho0=diag(diag(rho0));
end

function grumble(constants,paramsENDOR,B,geff,spinOps,spinSys,opt,HF_zz,HF_zy,HF_zx,NQI_zz)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isnumeric(B)
    error('B must be numeric.');
end
if ~isnumeric(geff)
    error('geff must be numeric.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
if ~isnumeric(HF_zz)
    error('HF_zz must be numeric.');
end
if ~isnumeric(HF_zy)
    error('HF_zy must be numeric.');
end
if ~isnumeric(HF_zx)
    error('HF_zx must be numeric.');
end
if ~isnumeric(NQI_zz)
    error('NQI_zz must be numeric.');
end
end

