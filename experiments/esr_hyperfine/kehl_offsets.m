% calculates the electron spin offsets for all spin manifolds in dependence
% of the hf values and nqc values
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% paramsENDOR: the Map containing the ENDOR parameters
% B: the main field
% geff: the effective g value
% HF_zz: the effective HF coupling value
% NQI_zz: the effective NQC value
% n: number of nucleus for multiple spin systems
%
% output parameters:
% offsets: double with offsets for all spin manifolds
%
% February 2024 A. Kehl (akehl@gwdg.de)
%



function offsets=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz)

    % Check consistency
    grumble(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);
    % get parameters out of the maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    v_L=paramsENDOR("v_L");
    I=spinSys("I");

    % loop over nuclei if not all in one spinSys
    for n=1: spinSys('N_SpinSys')
        % spin Hamiltonian
        H_EZ=B*geff*constants("MU_B")/constants("H")*Sz;
        H_NZ=zeros(size(Sz));
        H_HF=zeros(size(Sz));
        H_NQI=zeros(size(Sz));

        if spinSys('N_SpinSys')>1
                H_NZ=H_NZ-v_L(n)*Iz{1};
                H_HF=H_HF+2*pi*HF_zz(n)*(Sz*Iz{1}) ;
                H_NQI=H_NQI+pi*NQI_zz(n)*(3*Iz{1}*Iz{1}-I(n)*(I(n)+1)*eye(size(Sz)));
        else
            for mm=1:Ni_ENDOR
                H_NZ=H_NZ-v_L(mm)*Iz{mm};
                H_HF=H_HF+2*pi*HF_zz(mm)*(Sz*Iz{mm}) ;
                H_NQI=H_NQI+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Sz)));
            end
        end

        H_S=H_EZ+H_NZ+H_HF+H_NQI;

        [V_EPR,E_EPR]=eig(H_S);
        E_EPR=real(diag(E_EPR));

        %initiallyze freq matrices
        freq_EPR=zeros(1,length(V_EPR)*length(V_EPR));
        trans_prob_EPR=zeros(1,length(V_EPR)*length(V_EPR));


        % Calculate transitions between all elements
        q=1;
        for x=1:length(V_EPR)
            for y=x+1:length(V_EPR)
                % Electron transition probability in the eigenbasis
                trans_prob_EPR(q)=abs(round((V_EPR(:,x))'*Sx*(V_EPR(:,y)),9))^2;

                %in Hz, Transition Frequency
                freq_EPR(q)=abs(E_EPR(x)-E_EPR(y));

                q=q+1;
            end
        end

        freq_EPR=freq_EPR(logical(trans_prob_EPR));

        offsets_tmp=zeros(1,size(freq_EPR,2));

        % calculate offsets
        if mod(size(freq_EPR,2),2)==0
            for mm=1:size(freq_EPR,2)/2
                offsets_tmp(mm)=-(freq_EPR(size(freq_EPR,2)+1-mm)-freq_EPR(mm))/2;
            end
            for mm=size(freq_EPR,2)/2+1 : size(freq_EPR,2)
                offsets_tmp(mm)=-offsets_tmp(size(freq_EPR,2)+1-mm);
            end
        else
            center=size(freq_EPR,2)/2+0.5;
            for mm=1:size(freq_EPR,2)
                offsets_tmp(mm)=freq_EPR(mm)-freq_EPR(center);
            end
        end
        r=size(freq_EPR,2);

        offsets ((1+(n-1)*r):n*r)=offsets_tmp/(2*pi);

    end


end

function grumble(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
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
if ~isnumeric(HF_zz)
    error('HF_zz must be numeric.');
end
if ~isnumeric(NQI_zz)
    error('NQI_zz must be numeric.');
end
end

