% adds nuclei to be considered only for the EPR calculation to the spinsystem
% is directly called for the set up
%
% input parameters:
% constants: the Map with the constants
% spinSys_in: the Map with the spin parameters
% EPR_Nuclei: the nuclei for the EPR calculation (e.g. '14N')
% A_ten_EPR: the diagonal values of the hf tensor in MHz
% A_ang_EPR: the Euler angles for the hf tensor (alpha, beta, gamma)
% Q_ten_EPR: the diagonal values of the nq tensor in MHz
% Q_ang_EPR: the Euler angles for the nq tensor (alpha, beta, gamma)
%
% output parameters:
% SpinSys_out: updated Map with spin parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)

function SpinSys_out=kehl_spin_epr(constants,spinSys_in,EPR_Nuclei,A_ten_EPR,A_ang_EPR,Q_ten_EPR,Q_ang_EPR)
    if nargin<6
        Q_ten_EPR=[];
    end
    if nargin<7
        Q_ang_EPR=[];
    end


    % Check consistency
    grumble(constants,spinSys_in,EPR_Nuclei,A_ten_EPR,A_ang_EPR,Q_ten_EPR,Q_ang_EPR);
    SpinSys_out=spinSys_in;

    SpinSys_out("EPR_Nuclei")=EPR_Nuclei;
    Ni_EPR=size(EPR_Nuclei,2);
    SpinSys_out("Ni_EPR")=Ni_EPR;

    I=zeros(Ni_EPR,1);
    g_N_EPR=zeros(Ni_EPR,1);
    for i=1:Ni_EPR
        [~,multiplicity]=spin(kehl_spin_label(EPR_Nuclei{i}));
        I(i)=(multiplicity-1)/2;
        g_N_EPR(i)=kehl_nuc_gamma(constants,EPR_Nuclei{i});
    end
    SpinSys_out("I_EPR")=I;
    SpinSys_out("g_N_EPR")=g_N_EPR;

    SpinSys_out("EPR_Nucs_used")=true;
    SpinSys_out("EPR_Q_used")=false;

    L_EPR=zeros(Ni_EPR,12);
    for i=1:Ni_EPR
        L_EPR(i,1:3)=A_ten_EPR(3*i-2:3*i)*10^6;
        L_EPR(i,4:6)=A_ang_EPR(3*i-2:3*i);
        if size(Q_ten_EPR,2)>0
            SpinSys_out("EPR_Q_used")=true;
            L_EPR(i,7:9)=Q_ten_EPR(3*i-2:3*i)*10^6;
            L_EPR(i,10:12)=Q_ang_EPR(3*i-2:3*i);
        end
    end
    [A_EPR,Q_EPR]=kehl_to_g_frame(L_EPR,Ni_EPR);

    SpinSys_out("A_EPR")=A_EPR;
    SpinSys_out("Q_EPR")=Q_EPR;
end

function grumble(constants,spinSys_in,EPR_Nuclei,A_ten_EPR,A_ang_EPR,Q_ten_EPR,Q_ang_EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinSys_in,'containers.Map')
    error('spinSys_in must be a containers.Map object.');
end
if (~iscell(EPR_Nuclei))&&(~ischar(EPR_Nuclei))&&(~isstring(EPR_Nuclei))
    error('EPR_Nuclei must be a cell array, a character string, or a string scalar.');
end
if ~isnumeric(A_ten_EPR)
    error('A_ten_EPR must be numeric.');
end
if ~isnumeric(A_ang_EPR)
    error('A_ang_EPR must be numeric.');
end
if ~isnumeric(Q_ten_EPR)
    error('Q_ten_EPR must be numeric.');
end
if ~isnumeric(Q_ang_EPR)
    error('Q_ang_EPR must be numeric.');
end
end

