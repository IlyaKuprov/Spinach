% performs the actual ENDOR calculation (no relaxation)
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% opt: the Map containing the optional paramters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_mims_calc(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR)

    % Check consistency
    grumble(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR);
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,opt,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
                if opt("powder")==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if opt("Bterm")==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if opt("Bterm")==true
                        HNQI=zeros(size(Hfree_p));
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,1)*Ix{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,2)*Iy{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,3)*Iz{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,1)*Ix{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,2)*Iy{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,3)*Iz{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,1)*Ix{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,2)*Iy{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,3)*Iz{mm};
                        Hfree_p=Hfree_p+HNQI;
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));

                    end

                end

                if spinSys("D_used")==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;

                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % mw pulses

            %Hfree_p +
            Hnonsel_p=oneE*Sx;

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/abs(off_1*Nint);
            else
                t9=1/abs(v_off_S*Nint);
            end


            % Calculate the propagators
            U1_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(1)));
            U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));
            if t(3)==t(1)
                U3_p=U1_p;
            else
                U3_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(3)));
            end

            U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));

            if t(4)==t(6)
                U6_p=U4_p;
            else

                U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));
            end

            if t(7)==t(1)
                U7_p=U1_p;
            else
                U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
            end
            U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)+t(3)/2));
            U9_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if opt("Bterm")==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                    end
                end


                if opt("Bterm")==false
                    U5=full(propagator(spinOps('spin_system'),sparse(HRF),t(5)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));
                    if t(3)==t(1)
                        U3=U1;
                    else
                        U3=U3_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(3)));
                    end
                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));

                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));

                    if t(4)==t(6)
                        U6=U4;
                    else
                        U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));
                    end

                    if t(7)==t(1)
                        U7=U1;
                    else
                        U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));
                    end
                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)+t(3)/2));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U5=kehl_rf_bterm(opt,expt,v_RF,Hfree,Iy,t(5),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;
                    U3=U3_p;
                    U4=U4_p;

                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                end

                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                rho=U3*rho*U3';
                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if opt("Bterm")==true
                    rho=diag(diag(rho));
                end

                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';


                value_Sy=0;
                value_Sy=value_Sy+(real(trace(rho*Sy)));
                for b=1:Nint*10
                   rho=U9*rho*U9';
                   value_Sy=value_Sy+(real(trace(rho*Sy)));
                end

                % endor_amp_tmp(a) = endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(offsets,2)));

                if spinSys('N_SpinSys')>1
                    if opt('powder')==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end
                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*s));

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function grumble(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end

