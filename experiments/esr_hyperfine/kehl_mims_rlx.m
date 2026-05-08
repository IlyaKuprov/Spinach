% performs the actual ENDOR calculation (with relaxation)
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


function endor_amp=kehl_mims_rlx(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR)

    % Check consistency
    grumble(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");

    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");

    t=expt("t");
    Nint=100;


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


        if abs(HF_zz(1))>1/spinSys("T2e")*0.1

            % loop over all offsets (spin manifolds)
            for i=1 : size(offsets,2)

                v_off_S=offsets(i);
                off_1=offsets(1);

                rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,opt,HF_zz,HF_zy,HF_zx,NQI_zz);

                %electron T1
                RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
                RT1n=zeros(size(RT1e));
                RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
                RT2n=zeros(size(RT2e));

                if N_spinSys==1
                    for mm=1:Ni_ENDOR

                       %nuclear T1
                       RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));
                       RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                       RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                    end
                else
                    RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                    RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                    RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
                end
                R=RT1e+RT1n+RT2e+RT2n;


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
                            Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                            Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
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

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t9=1/(off_1*Nint);
                else
                    t9=1/(v_off_S*Nint);
                end

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
                            Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                            HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                        else
                            for mm=1:Ni_ENDOR
                                Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                                HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                            end
                        end
                    end

                    Hfree=Hfree_p+Hcorr;

                    Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                    Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    if opt("Bterm")==false

                        U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));
                    else
                        U5=kehl_rf_bterm_rlx(opt,expt,v_RF,Hfree_p,Iy,t(5),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                        U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));

                    end


                    % Evolve the densitymatrix
                    rho=kehl_mat_to_lbra(rho0)';

                    rho=U1*rho;
                    rho=U2*rho;
                    rho=U3*rho;
                    rho=U4*rho;
                    rho=U5*rho;
                    rho=U6*rho;
                    rho=U7*rho;
                    rho=U8*rho;


                    value_Sy=0;
                    for b=1
                       rho=U9*rho;
                       rho_f=kehl_lket_to_mat(rho);
                       value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                    end
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

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

