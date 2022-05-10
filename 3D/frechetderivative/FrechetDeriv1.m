function Jac = FrechetDeriv1(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,options,FwData)
%% Created on 19 November 2019
% Modified on 8th feb 2020.
% Now includes support for different phase functions


%
%% Evaluating the Frechet Derivative
Jac=[];
NofNodes = size(mesh.nodes,1);
NofElem = size(mesh.tri,1);
deflag=0;
if (isa(AssemblyFile,'function_handle'))
    disp('Obtain FEM matrices...')
    [Ks,Km,Kb,Iv,Jv] = AssemblyFile(mesh);
else
    load(AssemblyFile,'Ks','Km','Kb','Iv','Jv');
end
%     [~,e_n]=boundedges_element(P,T);    %boundary edges
Axx=AdjointFwd(Dimension,N,FwData.cgradphix,FwData.cphix,FwData.cgradbx,FwData.cbx,mesh,Ks,Km,Kb,Iv,Jv);
% Axx = FwData.sparsex';
NofDet=length(Det.elem);
NofSources = size(FwData.phix,2);
% onemat = ones(1,1,NofElem);
if(SimType)
    muax=zeros(1,1,NofElem);
    %set total muax, muam
    muax(1,1,:)=mesh.opt.muaxi(1,1,:)+mesh.opt.muaxf(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed);
    muam=zeros(1,1,NofElem); %absorption coefficient at emission wavelength
    muam(1,1,:)=mesh.opt.muaimult*mesh.opt.muaxi(1,1,:)+mesh.opt.muafmult*mesh.opt.muaxf(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed);
    %     Amm = FwData.sparsem';
    Amm=AdjointFwd(Dimension,N,FwData.cgradphim,FwData.cphim,FwData.cgradbm,FwData.cbm,mesh,AssemblyFile);
    if (JacType)
        Del=zeros(NofNodes,NofDet);
        for i=1:NofDet
            S=zeros(NofNodes,NofDet);
            S(Det.nodes(i))=1;
            Del(:,i)=S;
        end
        switch(N)
            case 3
                Del=[Del, zeros(NofNodes,NofDet); zeros(NofNodes,NofDet),Del];
            case 5
                Del=[Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),Del,zeros(NofNodes,NofDet);zeros(Nodes,NofDet),zeros(NofNodes,NofDet),Del];
            case 7
                Del=[Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del,zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),Del];
        end
    else
        Del=GetPartCurr(Dimension,N,FwData.cjm,mesh.nodes,mesh.tri,Det)';
    end
    MBetaT=GetAssembledMat(Dimension,N,Iv,Jv,NofNodes,permute(conj(FwData.mbeta),[2 1 3]),Km);
    
    switch options.itersolver
        case 'none'
            PsiMM=Amm\Del;
            PsiXM=Axx\(MBetaT*PsiMM);
        case 'gmres'
            parfor d=1:NofDet
                PsiMM(:,d) = gmres(Amm,Del(:,d),options.restart,options.tol, options.maxit);
                PsiXM(:,d) = gmres(Axx,MbetaT*PsiMM(:,d),options.restart,options.tol, options.maxit);
            end
        case 'bicgstab'
            parfor d=1:NofDet
                PsiMM(:,d) = bicgstab(Amm,Del(:,d),options.tol, options.maxit);
                PsiXM(:,d) = bicgstab(Axx,MbetaT*PsiMM(:,d),options.restart,options.tol, options.maxit);
            end
        case 'pcg'
            parfor d=1:NofDet
                PsiMM(:,d) = pcg(Amm,Del(:,d),options.tol, options.maxit);
                PsiXM(:,d) = pcg(Axx,MbetaT*PsiMM(:,d),options.restart,options.tol, options.maxit);
            end
    end
     
    %         PsiMM=(AM')\Del;
    %     load(AssemblyFile,'Iv','Jv','Km');
    %         PsiXM=(AX')\(MBetaT*PsiMM);
    
    clear Del JmatM JmatX
else
    muax=zeros(1,1,NofElem);
    %set total muax, muam
    muax(1,1,:)=mesh.opt.muaxi(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,NofElem);
    if (JacType)
        Del=zeros(NofNodes,NofDet);
        for i=1:NofDet
            S=zeros(NofNodes,1);
            S(Det.nodes(i))=1;
            Del(:,i)=S;
        end
        switch(N)
            case 3
                Del=[Del, zeros(NofNodes,NofDet); zeros(NofNodes,NofDet),Del];
            case 5
                Del=[Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),Del,zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),Del];
            case 7
                Del=[Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(Nodes,NofDet),Del,zeros(NofNodes,NofDet),zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),Del,zeros(NofNodes,NofDet);zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),zeros(NofNodes,NofDet),Del];
        end
    else
        Del=GetPartCurr(Dimension,N,FwData.cjx,mesh.nodes,mesh.tri,Det)';
    end
    switch options.itersolver
        case 'none'
            PsiXX=Axx\Del;
            
        case 'gmres'
            parfor d=1:NofDet
                PsiXX(:,d) = gmres(Axx,Del(:,d),options.restart,options.tol, options.maxit);
                
            end
        case 'bicgstab'
            parfor d=1:NofDet
                PsiXX(:,d) = bicgstab(Axx,Del(:,d),options.tol, options.maxit);
                
            end
        case 'pcg'
            parfor d=1:NofDet
                PsiXX(:,d) = pcg(Axx,Del(:,d),options.tol, options.maxit);
                
            end
    end
    
    clear Del JmatX
end
switch(options.pfuncflag)
    case 0
        bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx);
    case 1
        bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.gamma);
    case 2
        bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.alpha, mesh.opt.gbx);
end
[dCgradphiX,dCphiX,dCgradbphiX,~]=dSP3Formulation(N,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Det);
% [dCgradphiX,dCphiX,dCgradbphiX,~]=dSP3FormulationN(muax,musx,g_x,e_n,detE,Elem,nfibre,AOD); % obtain analytical derivative of coeffecient matrices

Cmat_zeros=zeros((N+1)/2,(N+1)/2,NofElem); % Zero matrix . More efficient to use this instead of using zeros() repeatedly in code.
dbound_greenX=Cmat_zeros;
dCbsX=Cmat_zeros;
dCbsX2=Cmat_zeros;
CgradbpX=Cmat_zeros;
%     if N>1
for i=1:size(mesh.edgeElem,1)
    dCbsX(:,:,mesh.edgeElem(i,1))=dCgradphiX(:,:,mesh.edgeElem(i,1))/FwData.cgradbx(:,:,mesh.edgeElem(i,1));
    dbound_greenX(:,:,mesh.edgeElem(i,1))=dCbsX(:,:,mesh.edgeElem(i,1))*FwData.cbx(:,:,mesh.edgeElem(i,1)); %Term from green's identity in FEM formulation
    
    dCbsX2(:,:,mesh.edgeElem(i,1))=(FwData.cgradphix(:,:,mesh.edgeElem(i,1))/FwData.cgradbx(:,:,mesh.edgeElem(i,1)))*dCgradbphiX(:,:,mesh.edgeElem(i,1));
    CgradbpX(:,:,mesh.edgeElem(i,1))=(dCbsX2(:,:,mesh.edgeElem(i,1))/FwData.cgradbx(:,:,mesh.edgeElem(i,1)))*FwData.cbx(:,:,mesh.edgeElem(i,1)); %Cgradbp in notes
end
clear dCbsX dCbsX2;

if (SimType)
    switch(options.pfuncflag)
        case 0
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm);
        case 1
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm,mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm, mesh.opt.alpha, mesh.opt.gbm);
    end
    [dCgradphiM,dCphiM,dCgradbphiM,~]=dSP3Formulation(N,mesh,muam,mesh.opt.musm,bn,mesh.opt.refrind_out,mesh.opt.refrind,Det,mesh.opt.muafmult);
    
    dbound_greenM=Cmat_zeros;
    dCbsM=Cmat_zeros;
    dCbsM2=Cmat_zeros;
    CgradbpM=Cmat_zeros;
    for i=1:size(mesh.edgeElem,1)
        dCbsM(:,:,mesh.edgeElem(i,1))=dCgradphiM(:,:,mesh.edgeElem(i,1))/FwData.cgradbm(:,:,mesh.edgeElem(i,1));
        dbound_greenM(:,:,mesh.edgeElem(i,1))=dCbsM(:,:,mesh.edgeElem(i,1))*FwData.cbm(:,:,mesh.edgeElem(i,1)); %Term from green's identity in FEM formulation
        
        dCbsM2(:,:,mesh.edgeElem(i,1))=(FwData.cgradphim(:,:,mesh.edgeElem(i,1))/FwData.cgradbm(:,:,mesh.edgeElem(i,1)))*dCgradbphiM(:,:,mesh.edgeElem(i,1));
        CgradbpM(:,:,mesh.edgeElem(i,1))=(dCbsM2(:,:,mesh.edgeElem(i,1))/FwData.cgradbm(:,:,mesh.edgeElem(i,1)))*FwData.cbm(:,:,mesh.edgeElem(i,1)); %Cgradbp in notes
    end
    clear dCbsM dCbsM2
    dmuaxf=ones(1,1,NofElem);
    db=(mesh.opt.Qf/(1-1i*Src(1).mfreq*mesh.opt.tau)).*dmuaxf;
    switch(N)
        case 1 
            sbeta = db; 
        case (3)
            dbeta=[db (-2/3)*db; (-2/3)*db (4/9)*db];
        case (5)
            dbeta=[db (-2/3)*db (8/15)*db
                (-2/3)*db (4/9)*db -(16/45)*db
                (8/15)*db -(16/45)*db (64/225*db)];
        case (7)
            dbeta=[db (-2/3)*db (8/15)*db -(16/35)*db
                (-2/3)*db (4/9)*db -(16/45)*db, (32/105)*db
                (8/15)*db -(16/45)*db (64/225*db) -(128/525)*db
                -(16/35)*db (32/105)*db -(128/525)*db (256/1225)*db];
    end
    
    clear db dmuaxf
end

D=Dimension+1;
% NofSources = size(Src,1);
LZ = zeros((N+1)/2*D,(N+1)/2*D,NofElem); %Zero matrix that can be used repeatedly
dLX=LZ;
PHIX=zeros(NofSources,(N+1)/2*D,NofElem);
TX = zeros((N+1)/2*D,NofSources,NofElem); %Zero matrix that can be used repeatedly

if(SimType)
    % Initialise matrices
    dLM=LZ;
    dLbeta=LZ;
    PHIM=PHIX;
    PSIXMT=zeros((N+1)/2*D,NofDet,NofElem);
    PSIMMT=PSIXMT;
    TM=TX;
    Tbeta=TM;
    % Elementwise storage of nodal quantities
    for e=1:NofElem
        %Not for N=1
        PHIX(:,(1:(N+1)/2:(N+1)/2*D),e)=FwData.phix(mesh.tri(e,:),:).';
        PHIM(:,(1:(N+1)/2:(N+1)/2*D),e)=FwData.phim(mesh.tri(e,:),:).';
        PSIXMT((1:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiXM(mesh.tri(e,:),1:NofDet));
        PSIMMT((1:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiMM(mesh.tri(e,:),1:NofDet));
        if N>1
            PHIX(:,(2:(N+1)/2:(N+1)/2*D),e)=FwData.phix(NofNodes+mesh.tri(e,:),:).';
            PHIM(:,(2:(N+1)/2:(N+1)/2*D),e)=FwData.phim(NofNodes+mesh.tri(e,:),:).';
            PSIXMT((2:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiXM(NofNodes+mesh.tri(e,:),1:NofDet));
            PSIMMT((2:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiMM(NofNodes+mesh.tri(e,:),1:NofDet));
        end
    end
else
    PSIXXT=zeros((N+1)/2*D,NofDet,NofElem);
    for e=1:NofElem
        PHIX(:,[1:(N+1)/2:(N+1)/2*D],e)=FwData.phix(mesh.tri(e,:),:).';
        PSIXXT([1:(N+1)/2:(N+1)/2*D],:,e)=conj(PsiXX(mesh.tri(e,:),1:NofDet));
        if N>1
            PHIX(:,[2:(N+1)/2:(N+1)/2*D],e)=FwData.phix(NofNodes+mesh.tri(e,:),:).';
            PSIXXT([2:(N+1)/2:(N+1)/2*D],:,e)=conj(PsiXX(NofNodes+mesh.tri(e,:),1:NofDet));
        end
    end
end
for s=1:NofSources
    %Vectorized implementation of the Jacobian computation ref. Fedele
    lim1=(N+1)/2; lim2=(N+1)/2-1;
    PHIX_s=PHIX(s,:,:);
    if(SimType)
        PHIM_s=PHIM(s,:,:);
    end
    for ii=1:D
        for jj=1:D
            dLX(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiX + repmat(Km(ii,jj,:),lim1,lim1).*dCphiX + repmat(Kb(ii,jj,:),lim1,lim1).*(dbound_greenX-CgradbpX);
            if (SimType)
                dLM(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiM + repmat(Km(ii,jj,:),lim1,lim1).*dCphiM + repmat(Kb(ii,jj,:),lim1,lim1).*(dbound_greenM - CgradbpM);
                dLbeta(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Km(ii,jj,:),lim1,lim1).*dbeta;
            end
        end
        if N==1
            TX(lim1*ii-lim2:lim1*ii,s,:)=sum(dLX(lim1*ii-lim2:lim1*ii,:,:).*PHIX_s,2);
            if(SimType)
                TM(lim1*ii-lim2:lim1*ii,s,:)=sum(dLM(lim1*ii-lim2:lim1*ii,:,:).*PHIM_s,2);
                Tbeta(lim1*ii-lim2:lim1*ii,s,:)=sum(dLbeta(lim1*ii-lim2:lim1*ii,:,:).*PHIX_s,2);
            end
        else if N==3
                TX(lim1*ii-lim2:lim1*ii,s,:)=sum(dLX(lim1*ii-lim2:lim1*ii,:,:).*[PHIX_s;PHIX_s],2);
                if(SimType)
                    TM(lim1*ii-lim2:lim1*ii,s,:)=sum(dLM(lim1*ii-lim2:lim1*ii,:,:).*[PHIM_s;PHIM_s],2);
                    Tbeta(lim1*ii-lim2:lim1*ii,s,:)=sum(dLbeta(lim1*ii-lim2:lim1*ii,:,:).*[PHIX_s;PHIX_s],2);
                end
            end
        end
    end
    
    if (SimType)
        delphi=-sum(PSIXMT.*repmat(TX(:,s,:),1,NofDet),1)+sum(PSIMMT.*repmat(Tbeta(:,s,:),1,NofDet),1)-sum(PSIMMT.*repmat(TM(:,s,:),1,NofDet),1);
    else
        delphi=-sum(PSIXXT.*repmat(TX(:,s,:),1,NofDet),1);
    end
    %         JacobianF(1:d,:)=delphi(1,1:d,:);
    
    JacobianJ(1:NofDet,:)=delphi(1,1:NofDet,:);
    Jac=[Jac;JacobianJ];
    %         JacF=[JacF;JacobianF];
end
