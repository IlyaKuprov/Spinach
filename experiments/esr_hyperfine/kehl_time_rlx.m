% performs the actual ENDOR calculation (with relaxation)
% CAUTION: this is a beta version and not fully tested
%
% input parameters:
% constants: the Map containing the constants
% spinOps: the Map containing the spin operators
% spinSys: the Map describing the spin system
% expt: the Map containing the experimental parameters
% parameters: structure containing simulation parameters
% paramsENDOR: the Map containing the ENDOR parameters
% EPR: the Map with results from EPR calculation
%
% output parameters:
% endor_amp: endor amplitude values
%
% April 2024 A. Kehl (akehl@gwdg.de)
%

% get values from Maps

function endor_amp=kehl_time_rlx(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)

    % Check consistency
    grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

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
    warning('This is a beta version and has not been fully tested.');


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        endor_amp_tmp=zeros(1,Npts_EN);

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


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mn=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mn},spinSys("T1n"));
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mn},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mn},spinSys("T2n"));
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
                nuc=round(i/(2*I(1)+1)) ;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(1)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(1,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(1)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mn=1:Ni_ENDOR
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{mn}+2*pi*HF_zz(mn)*(Sz*Iz{mn})+2*pi*HF_zy(mn)*(Sz*Iy{mn})+2*pi*HF_zx(mn)*(Sz*Ix{mn})+2*pi*v_L(mn)*CS_zz(mn)*Iz{mn};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{mn}*Iz{mn};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{mn}*Iz{mn}-I(mn)*(I(mn)+1)*eye(size(Hfree_p)));
                    end
                end

                if spinSys("D_used")==true
                    for mn=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mn-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mn};
                        HD=HD+Ix{1}*Iy{mn};
                        HD=HD+Ix{1}*Iz{mn};
                        HD=HD+Iy{1}*Ix{mn};
                        HD=HD+Iy{1}*Iy{mn};
                        HD=HD+Iy{1}*Iz{mn};
                        HD=HD+Iz{1}*Ix{mn};
                        HD=HD+Iz{1}*Iy{mn};
                        HD=HD+Iz{1}*Iz{mn};

                        HDip=dipC*(3*Iz{1}*Iz{mn}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t11=1/(off_1*Nint);
            else
                t11=1/(v_off_S*Nint);
            end

            % Radiofrequency
            v_RF=v_L(1);

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mn=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mn};

                            HRF=HRF+2*pi*v_RF*Iz{mn}+oneN*Iy{mn};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));



                if parameters.Bterm==false

                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(7)));
                    end

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));
                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));

                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                    end
                    U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                else
                    U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(7),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                    end

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(3)));

                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(8)));
                    end

                    U10=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                end

                % Evolve the densitymatrix
                rho=kehl_mat_to_lbra(rho0)';

                rho=U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho_t=rho;


            % loop over different separation of the RF pulses (x-axis)
            for a=1:Npts_EN
                t6=t(6)+step_EN*(a-1);
                U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t6));

                rho=rho_t;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;
                rho=U9*rho;
                rho=U9*rho;
                rho=U10*rho;


                value_Sy=0;
                for b=1
                   rho_f=kehl_lket_to_mat(rho);
                   value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                end

                if abs(HF_zz(1))>1/spinSys("T2e")*0.1
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function grumble(constants,spinOps,spinSys,expt,parameters,paramsENDOR,EPR)
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
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(EPR,'containers.Map')
    error('EPR must be a containers.Map object.');
end
end

