global Dimension N pfuncflag evalFluence probType JacType SimType AssemblyFile fcs_flag
simulation ='F';
SettingsFile = 'PhantomD1_S16';
ReconFile ='SettingsD1_S16';

norm_data.flag =0;
norm_data.ndet = 39;
for model =1
    % model =2;
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');    
     switch(model)
        case 1
            func = @FwdSolver;
            if norm_data.flag==1
                func_deriv = @FrechetDeriv1_mx;
            else
            func_deriv = @FrechetDeriv1;
            end
            fcs_flag=0;
            fprintf('Model: SPN-MFP \n');
        case 2
            func = @FwdSolver_fcs_5;
             if norm_data.flag==1
                 func_deriv = @FrechetDeriv1_fcs_mx_4;
             else
            func_deriv = @FrechetDeriv1_fcs_3;
             end
            fcs_flag=1;
            fprintf('Model: SPN-FCS \n');
        case 3
            func = @FwdSolver_fcs_de_4;
            
            if norm_data.flag==1
            func_deriv = @FrechetDeriv1_de_mx;
            else
                func_deriv = @FrechetDeriv1_de_3;
            end
            fcs_flag=1;
            fprintf('Model: SPN-DE \n');
    end


    
    JacType =0;
    
    snr = 100;
    fprintf('Generating measurement data with SNR = %d \n',snr)
    norm_flag=1; % flag for including normalisation
    if fcs_flag==1
        [mesh_orig,Md_n,phi_u]=Datagen(simulation,SettingsFile,func,snr,fcs_flag,norm_data.flag);
    else
        [mesh_orig,Md_n]=Datagen(simulation,SettingsFile,func,snr,fcs_flag,norm_data.flag);
    end
    
    
    %%
    fprintf('Setting up the inverse problem... \n');
    
    run(ReconFile) %Contains the Problem setting
    
    NofSources = size(Source,2);
    for s=1:NofSources
        if Source(s).AOI ~=90
            Source(s).AOI =pi - asin(Source(s).nfibre*sin(Source(s).AOI*pi/180)/Background_property(1,end));
            Source(s).AOI=Source(s).AOI*180/pi;
        end
        Source(s).mfreq = 2*pi*Source(s).mfreq*1e6;
    end
    
    if Dimension == 2 % Choose '2' for 2D and '3' for 3D
        if isempty(shape)
            [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector);
        else
            [mesh,Src,Det]=MeshGen2D_new(saveMesh,Domain,h,Source,Detector,shape,centre,r);
        end
        
    else
        % Not yet set up for 3D
        if isempty(shape)
            [mesh,Src,Det]=MeshGen3D(Domain,h,Src,Det);
        else
            [mesh,Src,Det]=MeshGen3D(Domain,h,Src,Det,shape,centre,r);
        end
    end
    
    if Dimension==3
        AssemblyFile=Assemble3D(mesh.nodes,mesh.tri);
    else if Dimension ==2
            AssemblyFile=Assemble2D_new(mesh);
        end
    end
    NofElem=size(mesh.tri,1);
    n0=external_ref_ind;
    
    if simulation =='F'
        mesh.opt.muaxi=Background_property(1)*ones(1,1,NofElem);
        mesh.opt.musx=Background_property(2)*ones(1,1,NofElem);
        mesh.opt.muaimult=Background_property(3);
        mesh.opt.musm=Background_property(4)*mesh.opt.musx;
        mesh.opt.gx=Background_property(5);
        mesh.opt.gm=Background_property(6);
        mesh.opt.refrind=Background_property(7);
        mesh.opt.speed=speed_of_light/mesh.opt.refrind;
        mesh.opt.refrind_out =n0;
        mesh.opt.gamma = gammaval;  % gamma parameter for modified HG
        mesh.opt.alpha = alphaval; % alpha parameter for TTHG
        mesh.opt.gbx = gb_x; % backscattering anisotropy fot TTHG
        mesh.opt.gbm = gb_m; % backscattering anisotropy fot TTHG
        %Data for fluorophore
        %     muaxf=zeros(1,1,NofElem);
        mesh.opt.muaxf=Fluorophore_property(1)*ones(1,1,NofElem);
        if ~isempty(shape)
            index_inhom2 = [0,cumsum(mesh.inhom.ind)];
            
            for ii=1:length(Inhom_property)
                %     mesh.opt.muaxf(1,1,tF(1+index_inhom2(ii)))= Inhom_property(1);
                mesh.opt.muaxf(1,1,mesh.inhom.elem(1+index_inhom2(ii):index_inhom2(ii+1)))=Inhom_property(ii);
            end
        end
        mesh.opt.muafmult=Fluorophore_property(2);
        mesh.opt.Qf=Fluorophore_property(4)*ones(1,1,NofElem);
        mesh.opt.tau=Fluorophore_property(5); % in ns
        
        %     Param={mesh.opt,Src,Det,nfibre};
        
        %     func = @FwdSolverF;
        SimType=1;
        %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxf);
        
    else
        %Data for Layer 1
        mesh.opt.muaxi=Background_property(1)*ones(1,1,NofElem);
        mesh.opt.musx=Background_property(2)*ones(1,1,NofElem);
        index_inhom2 = [0,cumsum(mesh.inhom.ind)];
        for ii=1:length(Inhom_property)
            mesh.opt.muaxi(1,1,mesh.inhom.elem(1+index_inhom2(ii):index_inhom2(ii+1)))=Inhom_property(ii);
        end
        mesh.opt.gx=Background_property(3);
        mesh.opt.refrind=Background_property(4);
        mesh.opt.speed=speed_of_light/mesh.opt.refrind;
        mesh.opt.refrind_out =n0;
        mesh.opt.gamma = gamma;  % gamma parameter for modified HG
        mesh.opt.alpha = alpha; % alpha parameter for TTHG
        mesh.opt.gbx = gbx; % backscattering anisotropy fot TTHG
        
        %Simulate for elastic scattering
        %     Param={mesh.opt,Src,Det};
        %     func = @FwdSolverE;
        SimType=0;
        %     Plot_mu(mesh.nodes,mesh.tri,mesh.hx,mesh.opt.muaxi);
    end
    if fcs_flag==1
        for s=1:NofSources
            temp =180-Src(s).AOI;
            if Src(s).Loc(1,1)==Domain(1)
                %         normal_in = '+x';
                %         Src(s).ang_in = 270-Src(s).AOI;
                Src(s).ang_in = 90+temp;
            elseif Src(s).Loc(1,1)==Domain(2)
                %         normal_in = '-x';
                Src(s).ang_in = 270+temp;%90-Src(s).AOI;
            elseif Src(s).Loc(1,2)==Domain(3)
                %         normal_in = '+z'
                Src(s).ang_in =temp; %180-Src(s).AOI ;
            elseif Src(s).Loc(1,2)==Domain(4)
                Src(s).ang_in=360-Src(s).AOI;
            end
            
            [Src(s).ray_elem,Src(s).seg_elem,Src(s).lseg_elem] = get_raytrace_2D(Domain,mesh.tri,mesh.nodes,Src(s).Loc(1,:),Src(s).elem,Src(s).ang_in*pi/180);
        end
    else
        mu_t_mean = mean(mesh.opt.muaxi + mesh.opt.muaxf + mesh.opt.musx*(1-mesh.opt.gx));
        mfp = 1/(mu_t_mean);
        % Adjust inocming angle and move sources 1mfp inside.
        for s=1:NofSources
            temp =180-Src(s).AOI;
            if Src(s).Loc(1,1)==Domain(1)
                %         normal_in = '+x';
                %         Src(s).ang_in = 270-Src(s).AOI;
                Src(s).ang_in = 90+temp;
                Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
                Src(s).ID = 'I';
                Src(s).AOI = 90;
                Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
                Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));

            elseif Src(s).Loc(1,1)==Domain(2)
                %         normal_in = '-x';
                Src(s).ang_in = 270+temp;%90-Src(s).AOI;
                Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
                Src(s).ID = 'I';
                Src(s).AOI = 90;
                Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
                Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
            elseif Src(s).Loc(1,2)==Domain(3)
                %         normal_in = '+z'
                Src(s).ang_in =temp; %180-Src(s).AOI ;
                Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
                Src(s).ID = 'I';
                Src(s).AOI = 90;
                Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
                Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
            elseif Src(s).Loc(1,2)==Domain(4)
                Src(s).ang_in=360-Src(s).AOI;
                Src(s).Loc = Src(s).Loc + [mfp*sind(Src(s).ang_in),mfp*cosd(Src(s).ang_in)];
                Src(s).ID = 'I';
                Src(s).AOI = 90;
                Src(s).elem = pointLocation(mesh.tri,Src(s).Loc(1),Src(s).Loc(2));
                Src(s).node = nearestNeighbor(mesh.tri,Src(s).Loc(:,1),Src(s).Loc(:,2));
            end
            
        end
    end
    
    %% Plot initial problem related data
    % Plot muaxf
    cs1 = 0.01; %cs2 = 0.75;
%     fig2 = Plot_mu2(mesh_orig.nodes,mesh_orig.tri,mesh_orig.hx,mesh_orig.hy,mesh_orig.opt.muaxf,4,cs1,cs2);
    fig2 = Plot_mu2(mesh_orig.nodes,mesh_orig.tri,mesh_orig.hx,mesh_orig.hy,mesh_orig.opt.muaxf,4,cs1);
    figure(fig2)
    hold on
    for s=1:NofSources
        scatter(Src(s).Loc(1,1),Src(s).Loc(1,2),'ro','filled')
    end
    NofDet = size(Det.Loc,1);
    for d=1:NofDet
        scatter(Det.Loc(d,1),Det.Loc(d,2),'kd','filled')
    end
    
    
    % probType ='recon';
    % [phi_u,Jd0,Fluence,FwData]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,evalFluence,probType);
    %% Obtain initial estimate
    muk_est=0.006;
    muk=muk_est*ones(1,1,NofElem);
    fprintf('Running inverse solver... \n')
    switch (model)
        case 1
%             fprintf('\n Basic SPN ....\n')
            tic
            nls_method =1; LS_type =1; update_type =2;
            L = get_incidence(mesh.tri); % set to L = zeros(NofElem) if no Laplacian is required
            LMData_spn=ResidualMinAlgo_optim(muk,Md_n,func,func_deriv,mesh,Src,Det,nls_method,LS_type,update_type,L,norm_data);
            LM_time = toc
            
            mu_LM=LMData_spn{1};
            nzs= length(find(LMData_spn{2})>0);
            temp_mu(1,1,:)=mu_LM(:,nzs);
%                     Plot_mu2(mesh.nodes,mesh.tri,h(1),h(2),temp_mu,1,cs1,cs2);
            Plot_mu2(mesh.nodes,mesh.tri,h(1),h(2),temp_mu,1,cs1);
            title('SPN-mfp')
            
        case 2
%             fprintf('\n SPN with FCS ....\n')
            tic
            nls_method =1; LS_type =1; update_type =2;
            L = get_incidence(mesh.tri); % set to L = zeros(NofElem) if no Laplacian is required
            LMData_fspn=ResidualMinAlgo_optim(muk,Md_n,func,func_deriv,mesh,Src,Det,nls_method,LS_type,update_type,L,norm_data);
            LM_time = toc
            
            mu_LM=LMData_fspn{1};
            nzs= length(find(LMData_fspn{2})>0);
            temp_mu(1,1,:)=mu_LM(:,nzs);
%                     Plot_mu2(mesh.nodes,mesh.tri,h(1),h(2),temp_mu,2,cs1,cs2);
            Plot_mu2(mesh.nodes,mesh.tri,h(1),h(2),temp_mu,2,cs1);
            title('SPN-fcs')
            
        case 3
%             fprintf('\n SPN with DE ....\n')
            tic
            nls_method =1; LS_type =1; update_type =2;
            L = get_incidence(mesh.tri); % set to L = zeros(NofElem) if no Laplacian is required
            LMData_dspn=ResidualMinAlgo_optim(muk,Md_n,func,func_deriv,mesh,Src,Det,nls_method,LS_type,update_type,L,norm_data);
            LM_time = toc
            
            mu_LM=LMData_dspn{1};
            nzs= length(find(LMData_dspn{2})>0);
            temp_mu(1,1,:)=mu_LM(:,nzs);
%             Plot_mu2(mesh.nodes,mesh.tri,h(1),temp_mu,3,cs1,cs2);
            Plot_mu2(mesh.nodes,mesh.tri,h(1),h(2),temp_mu,3,cs1);
            title('SPN-de')
            
    end
end