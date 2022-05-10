simulation ='F';
SettingsFile = 'PhantomD1_S16';
func = @FwdSolver;
fcs_flag=0;
run(SettingsFile) %Contains the Problem setting

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
probType ='recon';
% [phi_u,Jd0,Fluence,FwData]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType);
[Jd0,FwData]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType);
JacType=0;
% Jac = FrechetDeriv1_fcs_3(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% [JacX,JacM] = FrechetDeriv1_fcs_mx(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% [JacX,JacM] = FrechetDeriv1_fcs_mx_4(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% [JacX,JacM] = FrechetDeriv1_mx(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% JacM = FrechetDeriv1_de_3(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
JacM = FrechetDeriv1(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);

disp('Adjoint Frechet derivative evaluated...')


pert=1e-6;
muaxf_old = mesh.opt.muaxf;
NofDet = size(Jd0{1}(:,2),1);
FD = zeros(NofSources*NofDet,size(mesh.tri,1));
disp('begin finite difference evaluation...')
probType = 'fwd';
for e =1:size(mesh.tri,1)
%     if sum(ismember(e,Src.ray_elem)||ismember(e,Det.elem),1)==0
    
    mesh.opt.muaxf = muaxf_old;
    mesh.opt.muaxf(1,1,e) = mesh.opt.muaxf(1,1,e) + pert;
%     [phi_u_pert,Jd_pert,~]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType);
    [Jd_pert]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType);
    FDM(1:NofDet,e) =(Jd_pert{1}(:,2) - Jd0{1}(:,2))/pert;
%     FDX(1:NofDet,e) =(Jd_pert{1}(:,1) - Jd0{1}(:,1))/pert;
%     FDU(:,e) = (phi_u_pert.nodal - phi_u.nodal)/pert;
    
%     FD(NofDet+1:2*NofDet,e) =(Jd_pert{2}(:,2) - Jd0{2}(:,2))/pert;
%     FD(2*NofDet+1:3*NofDet,e) =(Jd_pert{3}(:,2) - Jd0{3}(:,2))/pert;
%     FD(3*NofDet+1:4*NofDet,e) =(Jd_pert{4}(:,2) - Jd0{4}(:,2))/pert;
%     end
end
%%
% JacX2 = JacX; JacM2 = JacM;
% JacX(:,Src.ray_elem)=0; JacM(:,Src.ray_elem)=0;
% JacX(:,Det.elem)=0; JacM(:,Det.elem)=0;
disp('Plotting...')
figure
subplot(1,3,1)
imagesc(abs(diag(Jd0{1}(:,2))\JacM - diag(Jd0{1}(:,1))\JacX))
colorbar
title('adjoint-scaled')
subplot(1,3,2)
imagesc(abs(diag(Jd0{1}(:,2))\FDM - diag(Jd0{1}(:,1))\FDX))
colorbar
title('fd-scaled')
subplot(1,3,3)
imagesc(abs(FDX)-abs(JacX))
colorbar
figure(gcf)

figure
subplot(1,3,1)
imagesc(abs(JacM))
colorbar
title('adjoint-m')
subplot(1,3,2)
imagesc(abs(FDM))
colorbar
title('fd-m')
subplot(1,3,3)
imagesc(abs(FDM)-abs(JacM))
colorbar
figure(gcf)

figure
subplot(1,3,1)
imagesc((abs(JacX)))
colorbar
title('adjoint-x')
subplot(1,3,2)
imagesc((abs(FDX)))
colorbar
title('fd-x')
subplot(1,3,3)
imagesc((abs(FDX))-(abs(JacX)))
colorbar
figure(gcf)