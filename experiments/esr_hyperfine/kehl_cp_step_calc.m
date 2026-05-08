% performs the actual CP step calculation
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
% rho: the final spin density matrix
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function rho=kehl_cp_step_calc(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR)

    % Check consistency
    grumble(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR);
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));

    % performs the actual calculation
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
    Npts_CP=expt("Npts_CP");

    if expt("Npts_CP")>1
        step_CP=expt("range_CP")/(expt("Npts_CP")-1);
    else
        step_CP=0;
    end

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

    offsets_sel=EPR("offsets");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)

        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;

        offsets=offsets_sel(j,:);

            % loop over nuclei if not all in one spinSys
        for n=1 : nucs

            % loop over all offsets (spin manifolds)
            for i=1 : size(offsets,2)

                v_off_S=offsets(i);
                off_1=offsets(1);

                if opt("powder")==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(opt,paramsENDOR,expt,spinSys,v_off_S,HF_zz,NQI_zz,n,i);
                end

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


                prep=expt("prep");
                sl=expt("SL");
                cp=expt("CP");

                Hfree_p=2*pi*v_off_S*Sz;
                if spinSys('N_SpinSys')>1
                    nuc=round(i/(2*I(1)+1)) ;
                    mm=1;
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*v_L(mm)*Iz{mm}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{mm})+2*pi*HF_zy(nuc)*(Sz*Iy{mm})+2*pi*HF_zx(nuc)*(Sz*Ix{mm});
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
                end

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t11=1/(off_1*Nint);
                else
                    t11=1/(v_off_S*Nint);
                end

                % loop over all cp freq (y-axis)
                for c=1:Npts_CP

                    if opt("powder")==true
                        v_CP=expt("start_CP")+step_CP*(c-1);
                    end

                    Hcorr=zeros(size(Hfree_p));
                    HSL=Hfree_p+sl*Sy;

                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_CP*Iz{m};
                    else

                        for mm=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                        end
                    end
                    HSL_c=HSL+Hcorr;
                    [V,D]=eig(HSL_c);
                    if ~issorted(diag(D))
                      [d,inds]=sort(diag(D));
                      D=diag(d);
                      V=V(:, inds);
                    end
                    %eigenvalues are order from negative large to positive large!
                    W=V^(-1);


                    if opt("Bterm")==false
                        if spinSys('N_SpinSys')>1
                            m=1;
                            HSL=HSL+2*pi*v_CP*Iz{m}+cp*Ix{m};
                        else
                            for mm=1:Ni_ENDOR
                                HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                            end
                        end
                    end



                    Hfree=Hfree_p+Hcorr;

                    Hprep=full(hilb2liouv(sparse(Hfree+prep*Sx),'comm'));

                    Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    HSL_t=W*HSL*W^(-1);
                    if opt("Bterm")==false
                        if opt("Relax_step")==0
                            U3=full(propagator(spinOps('spin_system'),1i*sparse(-1i*full(hilb2liouv(sparse(HSL_t),'comm'))),t(3)));
                            U1=full(propagator(spinOps('spin_system'),1i*sparse(-1i*Hprep),t(1)));
                            U2=full(propagator(spinOps('spin_system'),1i*sparse(-1i*Hfree),t(2)));
                        else
                            U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HSL_t),'comm'))),t(3)));
                            U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                            U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                        end
                    else
                        U3=kehl_cp_bterm_rlx(opt,expt,v_CP,HSL_t,Ix,t(3),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));
                        U1=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));
                    end


                    % Evolve the densitymatrix

                    rho=kehl_mat_to_lbra(rho0)';

                    rho=U1*rho;
                    rho=U2*rho;
                    rho_t=kehl_lket_to_mat(rho);
                    rho_t=W*rho_t*W^(-1);
                    rho=kehl_mat_to_lbra(rho_t)';
                    rho=U3*rho;

                    rho=kehl_lket_to_mat(rho);

                end
            end
        end
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

