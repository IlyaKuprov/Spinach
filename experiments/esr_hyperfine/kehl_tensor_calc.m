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

% redefine, one spin system is not possible

function endor_amp=kehl_tensor_calc(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR)

    % Check consistency
    grumble(constants,spinOps,spinSys,expt,opt,paramsENDOR,EPR);
    spinSys('N_SpinSys')=spinSys("Ni_ENDOR");
    spinOps=kehl_spin_ops(spinSys("S"),spinSys("I"),spinSys("Ni_ENDOR"),spinSys("Ni_ENDOR"));

    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("Ni_ENDOR");
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
        const_R=1/size(Sz,2)*constants('GE')*B/(2*pi*constants('K_B')*opt("T"));

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

        % loop over nuclei
        for i=1 : Ni_ENDOR
            I=spinSys('I');
            mI=-I(i):1:I(i);

            for jj=mI

                v_off_S=jj*HF_zz(i);
                off_1=-I(i)*HF_zz(i);

                [rho0]=2*Sz*Iz{i};

                start_EN=paramsENDOR("start_EN");
                step_EN=paramsENDOR("step_EN");

                oneE=expt("oneE");
                oneN=expt("oneN");

                Hfree_p=2*pi*v_off_S*Sz;
                nuc=i
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(nuc)*Iz{1}+2*pi*v_L(nuc)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if opt("Bterm")==true
                    m=1;
                    Hfree_p=Hfree_p+NQI(m,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(m,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(m,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(m,3,3)*Iz{1}*Iz{1};
                else
                    m=1;
                    Hfree_p=Hfree_p+pi*NQI_zz(m)*(3*Iz{m}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                if spinSys("D_used")==true && i==1
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

                % mw pulses
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t2=1/abs(off_1*Nint);
                else
                    t2=1/abs(v_off_S*Nint);
                end



                U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t2));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if opt("Bterm")==false
                   HRF=HRF+2*pi*v_RF*Iz{1}+oneN*Iy{1};
                end

                Hcorr=Hcorr+2*pi*v_RF*Iz{1};


                if opt("Bterm")==false
                    U1=full(propagator(spinOps('spin_system'),sparse(HRF),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t2));
                else
                    Hfree=Hfree_p;
                    U1=kehl_rf_bterm(opt,expt,v_RF,Hfree,Iy,t(1),Ni_ENDOR,N_spinSys,spinOps('spin_system'));
                    U2=U2_p;
                end


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';

                if opt("Bterm")==true
                    rho=diag(diag(rho));
                end


                value_Sy=0;
                value_Sy=value_Sy+(real(trace(rho*Sy)));
                for b=1:Nint*10
                   rho=U2*rho*U2';
                   value_Sy=value_Sy+(real(trace(rho*Sz*Iz{1})));
                end

                if opt('temp_eff')==true
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)))*const_R;
                else
                    endor_amp_tmp(a)=endor_amp_tmp(a)+abs(value_Sy*S/(Nint*size(mI,2)));

                end

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

