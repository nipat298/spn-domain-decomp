function [Jd,varargout]=FwdSolver(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,pfuncflag,evalFluence,probType)
%%
%Last modified on 9th April 2019
% Modified on 8th feb 2020.
% Now includes support for different phase functions

% Forward solver based on the SPN approximation to the radiative transfer
% equation. can support N = 1,3,5,7. Not tested for N = 1,5,7 for this
% version yet. Earlier versions offer full support.
% Basic SPN (BSPN) version
%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%
% Dimesnions - Choose 2 for 2D or 3 for 3D
% N  - order of the SPN approximation
% mesh - structure type variable containing meshing information
% Src - structure type variable containing the source information
% Det - structure type variable containing the detector information
% AssemblyFile - name of the mat file containing the mass,stiffness
% matrices and indexing vectors
% evalFluence - Flag set to '0' to suppress fluence output, else set to '1' to
% probType - 'fwd' for data generation, forward solve only,
%            'recon' for generating data in a reconstruction routine. saves
%            some variables that are needed in the reconstruction routine
%            as well as evaluation of adjoint sensitivity

%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%
% Jd - Exiting partial current
% if evalFluence ==1
% varargout{1} = Fluence
% varargout{2} = FwData (if probType = 'recon') - structure containing
% information required for sensitivity and reconstruction routine
% if evalFluence == 0 and probType = 'recon'
% varargout{1} = FwData
%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com

%%
% DT - node connectivity list
% P - coordinates of all nodes
% e_n - [list of edge elements, type of edge element] % details in function
% boundedges_elem
Elem=size(mesh.tri,1);
Nodes=size(mesh.nodes,1);
deflag=0;
if Dimension ==2
    T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3)];
else
    T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3),mesh.tri(:,4)];
end




NofSources = size(Src,2);
if (SimType)
    % Fluorescence
    
    %Initialise matrices
    muax=zeros(1,1,Elem);
    %set total muax, muam
    muax(1,1,:)=mesh.opt.muaxi(1,1,:)+mesh.opt.muaxf(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,Elem);
    
    switch(pfuncflag)
        case 0
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx);
        case 1
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx,mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx,mesh.opt.alpha,mesh.opt.gbx);
            
    end
    [CgradphiX,CphiX,CgradbphiX,CbphiX,CJX]=SP3Formulation(N,Elem,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for excitation
    
    %Solve the forward problem to compute PhiX and PhiM
    AX=Spnfwd(Dimension,N,CgradphiX,CphiX,CgradbphiX,CbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,AssemblyFile);
    for s=1:NofSources
        SX(:,s)=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,mesh.edgeType,Src(s));
    end
    %     SX = ones(2*Nodes,1);
    PhiX=AX\SX;
    
    % Emission field
    muam=zeros(1,1,Elem); %absorption coefficient at emission wavelength
    muam(1,1,:)=mesh.opt.muaimult*mesh.opt.muaxi(1,1,:)+mesh.opt.muafmult*mesh.opt.muaxf(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,Elem);
    
    switch(pfuncflag)
        case 0
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gm);
        case 1
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gm,mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gm, mesh.opt.alpha, mesh.opt.gbm);
    end
    [CgradphiM,CphiM,CgradbphiM,CbphiM,CJM]=SP3Formulation(N,Elem,mesh,muam,mesh.opt.musm,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for emission
    mesh.opt.beta=(mesh.opt.Qf/(1-1i*Src(1).mfreq*mesh.opt.tau)).*mesh.opt.muaxf(1,1,:); %(muax0(1,1,:) + muaxf(1,1,:)); %Note this is muaxf and not muax as in earlier versions
    switch(N)
        case (3)
            Mbeta=[mesh.opt.beta,-(2/3)*mesh.opt.beta; (-2/3)*mesh.opt.beta 4/9*mesh.opt.beta];
        case (5)
            Mbeta=[mesh.opt.beta,-(2/3)*mesh.opt.beta, (8/15)*mesh.opt.beta
                (-2/3)*mesh.opt.beta, 4/9*mesh.opt.beta, -(16/45)*mesh.opt.beta
                (8/15)*mesh.opt.beta, -(16/45)*mesh.opt.beta, (64/225)*mesh.opt.beta];
        case (7)
            Mbeta=[mesh.opt.beta,-(2/3)*mesh.opt.beta, (8/15)*mesh.opt.beta, -(16/35)*mesh.opt.beta
                (-2/3)*mesh.opt.beta, 4/9*mesh.opt.beta, -(16/45)*mesh.opt.beta, (32/105)*mesh.opt.beta
                (8/15)*mesh.opt.beta, -(16/45)*mesh.opt.beta, (64/225)*mesh.opt.beta, -(128/525)*mesh.opt.beta
                -(16/35)*mesh.opt.beta, (32/105)*mesh.opt.beta, -(128/525)*mesh.opt.beta, (256/1225)*mesh.opt.beta ];
    end
    AM=Spnfwd(Dimension,N,CgradphiM,CphiM,CgradbphiM,CbphiM,mesh.nodes,mesh.tri,mesh.edgeElem,AssemblyFile);
    % Obtain emission source
    load(AssemblyFile,'Iv','Jv','Km');
    SM = GetAssembledMat(Dimension,N,Iv,Jv,Nodes,Mbeta,Km)*PhiX;
    
    PhiM=AM\SM;
    
    if (evalFluence)
        switch(N)
            case 3
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:end,:);
                FM = PhiM(1:Nodes,:) - (2/3)*PhiM(Nodes+1:end,:);
            case 5
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:end,:);
                FM = PhiM(1:Nodes,:) - (2/3)*PhiM(Nodes+1:2*Nodes,:) + (8/15)*PhiM(2*Nodes+1:end,:);
            case 7
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:3*Nodes,:) - (16/35)*PhiX(3*Nodes+1:end,:);
                FM = PhiM(1:Nodes,:) - (2/3)*PhiM(Nodes+1:2*Nodes,:) + (8/15)*PhiM(2*Nodes+1:3*Nodes,:) - (16/35)*PhiM(3*Nodes+1:end,:);
        end
    end
    
    % Note: Partial Current computations are only for Surface detectors.
    % For an internal detector, the default measurement is fluence;
    if (Det.ID=='S')
        JmatX=GetPartCurr(Dimension,N,CJX,mesh.nodes,T,mesh.edges,mesh.edgeElem,mesh.edgeType,Det);
        JX=JmatX*PhiX;
        JmatM=GetPartCurr(Dimension,N,CJM,mesh.nodes,T,mesh.edges,mesh.edgeElem,mesh.edgeType,Det);
        JM=JmatM*PhiM;
    else
        JX = GetDetFluence(Det,PhiX);
        JM = GetDetFluence(Det,PhiM);
    end
    
    
    %     if Dimension ==2
    %         JmatX=GetExitCurrent2D(N,CJX,mesh.nodes,T,mesh.edgeElem,mesh.edgeType,Det);
    %         JX=JmatX*PhiX;
    %         JmatM=GetExitCurrent2D(N,CJM,mesh.nodes,T,mesh.edgeElem,mesh.edgeType,Det);
    %         JM=JmatM*PhiM;
    %     else
    %         JmatX=GetExitCurrent3DN(CJX,detE,mesh.nodes,T,mesh.edgeElem,DetLoc);
    %         JX=JmatX*PhiX;
    %         JmatM=GetExitCurrent3DN(CJM,detE,mesh.nodes,T,mesh.edgeElem,DetLoc);
    %         JM=JmatM*PhiM;
    %     end
    
    %     Jd=[JX(:),JM(:)];
    
    Jd = cell(NofSources,1);
    for s=1:NofSources
        Jd{s} = [JX(:,s),JM(:,s)];
    end
    if (evalFluence)
        varargout{1}=[FX;FM];
    end
    
    if strcmp(probType,'recon')
        FwData.cgradphix= CgradphiX;
        FwData.cphix= CphiX;
        FwData.cgradbx= CgradbphiX;
        FwData.cbx= CbphiX;
        FwData.cjx = CJX;
        FwData.sparsex = AX;
        FwData.phix=PhiX;
        
        FwData.cgradphim= CgradphiM;
        FwData.cphim= CphiM;
        FwData.cgradbm= CgradbphiM;
        FwData.cbm= CbphiM;
        FwData.cjm = CJM;
        FwData.mbeta=Mbeta;
        FwData.sparsem = AM;
        FwData.phim=PhiM;
        if (evalFluence)
            varargout{2} = FwData;
        else
            varargout{1} = FwData;
        end
    end
    
    
    
else
    % Elastic scattering case
    
    muax=zeros(1,1,Elem);
    muax(1,1,:)=mesh.opt.muaxi(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,Elem);
    switch(pfuncflag)
        case 0
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx);
        case 1
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx,mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(pfuncflag,deflag,mesh.opt.gx,mesh.opt.alpha, mesh.opt.hx);
    end
    [CgradphiX,CphiX,CgradbphiX,CbphiX,CJX]=SP3Formulation(N,Elem,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for excitation
    
    
    %Solve the forward problem to compute PhiX and PhiM
    AX=Spnfwd(Dimension,N,CgradphiX,CphiX,CgradbphiX,CbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,AssemblyFile);
    
    for s=1:NofSources
        SX(:,s)=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,mesh.edgeType,Src(s));
    end
    
    PhiX=AX\SX;
    
    JmatX=GetPartCurr(Dimension,N,CJX,mesh.nodes,T,mesh.edges,mesh.edgeElem,mesh.edgeType,Det);
    JX =JmatX*PhiX;
     Jd = cell(NofSources,1);
    for s=1:NofSources
    Jd{s} = JX(:,s);
    end
    
    if (evalFluence)
        switch(N)
            case 3
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:end,:);
            case 5
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:end,:);
            case 7
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:3*Nodes) - (16/35)*PhiX(3*nodes+1:end,:);
        end
    end
    if (evalFluence)
        varargout{1}=FX;
    end
    
    if strcmp(probType,'recon')
        FwData.cgradphix= CgradphiX;
        FwData.cphix= CphiX;
        FwData.cgradbx= CgradbphiX;
        FwData.cbx= CbphiX;
        FwData.cjx = CJX;
        FwData.sparsex = AX;
        FwData.phix=PhiX;
        if (evalFluence)
            varargout{2} = FwData;
        else
            varargout{1} = FwData;
        end
    end
    
end


% %% Evaluating the Frechet Derivative
% if (mode=='R')||(mode=='H')
%
%     if isempty(varargin)
%         JacType=0; % 0 -Jacobian w.rt partial current, 1 - Jacobian w.r.t fluence
%     else
%         JacType=varargin{1};
%     end
%     JacF=[]; JacJ=[];
%     load(AssemblyFile)
%     %     [~,e_n]=boundedges_element(P,T);    %boundary edges
%     Axx=AdjointFwdN(CgradphiX,CphiX,CgradbphiX,CbphiX,e_n,Nodes,Elem,AssemblyFile);
%     NofDet=length(DNodes);
%     if(SimType)
%         Amm=AdjointFwdN(CgradphiM,CphiM,CgradbphiM,CbphiM,e_n,Nodes,Elem,AssemblyFile);
%         if (JacType)
%             Del=zeros(Nodes,NofDet);
%             for i=1:NofDet
%                 S=zeros(Nodes,NofDet);
%                 S(DNodes(i))=1;
%                 Del(:,i)=S;
%             end
%             switch(N)
%                 case 3
%                     Del=[Del, zeros(Nodes,NofDet); zeros(Nodes,NofDet),Del];
%                 case 5
%                     Del=[Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del];
%                 case 7
%                     Del=[Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del];
%             end
%         else
%             Del=JmatM';
%         end
%         PsiMM=Amm\Del;
% %         PsiMM=(AM')\Del;
%         MBetaT=GetAssembledMatN(AssemblyFile,Nodes,permute(conj(Mbeta),[2 1 3]),Km);
% %         PsiXM=(AX')\(MBetaT*PsiMM);
%         PsiXM=Axx\(MBetaT*PsiMM);
%         clear Del JmatM JmatX
%     else
%         if (JacType)
%             Del=zeros(Nodes,NofDet);
%             for i=1:NofDet
%                 S=zeros(Nodes,1);
%                 S(DNodes(i))=1;
%                 Del(:,i)=S;
%             end
%             switch(N)
%                 case 3
%                     Del=[Del, zeros(Nodes,NofDet); zeros(Nodes,NofDet),Del];
%                 case 5
%                     Del=[Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del];
%                 case 7
%                     Del=[Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet),zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del,zeros(Nodes,NofDet);zeros(Nodes,NofDet),zeros(Nodes,NofDet),zeros(Nodes,NofDet),Del];
%             end
%         else
%             Del=JmatX';
%         end
%         PsiXX=(AX')\Del;
%         clear Del JmatX
%     end
%
%     [dCgradphiX,dCphiX,dCgradbphiX,~]=dSP3FormulationN(muax,musx,g_x,e_n,detE,Elem,nfibre,AOD); % obtain analytical derivative of coeffecient matrices
%
%     Cmat_zeros=zeros((N+1)/2,(N+1)/2,Elem); % Zero matrix . More efficient to use this instead of using zeros() repeatedly in code.
%     dbound_greenX=Cmat_zeros;
%     dCbsX=Cmat_zeros;
%     dCbsX2=Cmat_zeros;
%     CgradbpX=Cmat_zeros;
%     %     if N>1
%     for i=1:size(e_n,1)
%         dCbsX(:,:,e_n(i,1))=dCgradphiX(:,:,e_n(i,1))/CgradbphiX(:,:,e_n(i,1));
%         dbound_greenX(:,:,e_n(i,1))=dCbsX(:,:,e_n(i,1))*CbphiX(:,:,e_n(i,1)); %Term from green's identity in FEM formulation
%
%         dCbsX2(:,:,e_n(i,1))=(CgradphiX(:,:,e_n(i,1))/CgradbphiX(:,:,e_n(i,1)))*dCgradbphiX(:,:,e_n(i,1));
%         CgradbpX(:,:,e_n(i,1))=(dCbsX2(:,:,e_n(i,1))/CgradbphiX(:,:,e_n(i,1)))*CbphiX(:,:,e_n(i,1)); %Cgradbp in notes
%     end
%     clear dCbsX dCbsX2;
%
%     if (SimType)
%         [dCgradphiM,dCphiM,dCgradbphiM,~]=dSP3FormulationN(muam,musm,g_m,e_n,detE,Elem,nfibre,AOD,muaf_mult);
%         dbound_greenM=Cmat_zeros;
%         dCbsM=Cmat_zeros;
%         dCbsM2=Cmat_zeros;
%         CgradbpM=Cmat_zeros;
%         for i=1:size(e_n,1)
%             dCbsM(:,:,e_n(i,1))=dCgradphiM(:,:,e_n(i,1))/CgradbphiM(:,:,e_n(i,1));
%             dbound_greenM(:,:,e_n(i,1))=dCbsM(:,:,e_n(i,1))*CbphiM(:,:,e_n(i,1));
%
%             dCbsM2(:,:,e_n(i,1))=(CgradphiM(:,:,e_n(i,1))/CgradbphiM(:,:,e_n(i,1)))*dCgradbphiM(:,:,e_n(i,1));
%             CgradbpM(:,:,e_n(i,1))=(dCbsM2(:,:,e_n(i,1))/CgradbphiM(:,:,e_n(i,1)))*CbphiM(:,:,e_n(i,1));
%         end
%         clear dCbsM dCbsM2
%         dmuaxf=ones(1,1,Elem);
%         db=(Q_f/(1-1i*w*tau)).*dmuaxf;
%         switch(N)
%             case (3)
%                 dbeta=[db (-2/3)*db; (-2/3)*db (4/9)*db];
%             case (5)
%                 dbeta=[db (-2/3)*db (8/15)*db
%                     (-2/3)*db (4/9)*db -(16/45)*db
%                     (8/15)*db -(16/45)*db (64/225*db)];
%             case (7)
%                 dbeta=[db (-2/3)*db (8/15)*db -(16/35)*db
%                     (-2/3)*db (4/9)*db -(16/45)*db, (32/105)*db
%                     (8/15)*db -(16/45)*db (64/225*db) -(128/525)*db
%                     -(16/35)*db (32/105)*db -(128/525)*db (256/1225)*db];
%         end
%
%         clear db dmuaxf
%     end
%
%     D=Dimension+1;
%     NofSources = size(Src,1);
%     LZ = zeros((N+1)/2*D,(N+1)/2*D,Elem); %Zero matrix that can be used repeatedly
%     dLX=LZ;
%     PHIX=zeros(NofSources,(N+1)/2*D,Elem);
%     TX = zeros((N+1)/2*D,NofSources,Elem); %Zero matrix that can be used repeatedly
%
%     if(SimType)
%         % Initialise matrices
%         dLM=LZ;
%         dLbeta=LZ;
%         PHIM=PHIX;
%         PSIXMT=zeros((N+1)/2*D,NofDet,Elem);
%         PSIMMT=PSIXMT;
%         TM=TX;
%         Tbeta=TM;
%         % Elementwise storage of nodal quantities
%         for e=1:Elem
%             %Not for N=1
%             PHIX(:,(1:(N+1)/2:(N+1)/2*D),e)=PhiX(T(e,:),:).';
%             PHIM(:,(1:(N+1)/2:(N+1)/2*D),e)=PhiM(T(e,:),:).';
%             PSIXMT((1:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiXM(T(e,:),1:NofDet));
%             PSIMMT((1:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiMM(T(e,:),1:NofDet));
%             if N>1
%                 PHIX(:,(2:(N+1)/2:(N+1)/2*D),e)=PhiX(Nodes+T(e,:),:).';
%                 PHIM(:,(2:(N+1)/2:(N+1)/2*D),e)=PhiM(Nodes+T(e,:),:).';
%                 PSIXMT((2:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiXM(Nodes+T(e,:),1:NofDet));
%                 PSIMMT((2:(N+1)/2:(N+1)/2*D),:,e)=conj(PsiMM(Nodes+T(e,:),1:NofDet));
%             end
%         end
%     else
%         PSIXXT=zeros((N+1)/2*D,NofDet,Elem);
%         for e=1:Elem
%             PHIX(:,[1:(N+1)/2:(N+1)/2*D],e)=PhiX(T(e,:),:).';
%             PSIXXT([1:(N+1)/2:(N+1)/2*D],:,e)=conj(PsiXX(T(e,:),1:NofDet));
%             if N>1
%                 PHIX(:,[2:(N+1)/2:(N+1)/2*D],e)=PhiX(Nodes+T(e,:),:).';
%                 PSIXXT([2:(N+1)/2:(N+1)/2*D],:,e)=conj(PsiXX(Nodes+T(e,:),1:NofDet));
%             end
%         end
%     end
%     for s=1:NofSources
%         %Vectorized implementation of the Jacobian computation ref. Fedele
%         lim1=(N+1)/2; lim2=(N+1)/2-1;
%         PHIX_s=PHIX(s,:,:);
%         if(SimType)
%             PHIM_s=PHIM(s,:,:);
%         end
%         for ii=1:D
%             for jj=1:D
%                 dLX(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiX + repmat(Km(ii,jj,:),lim1,lim1).*dCphiX + repmat(Kb(ii,jj,:),lim1,lim1).*(dbound_greenX -CgradbpX);
%                 if (SimType)
%                     dLM(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiM + repmat(Km(ii,jj,:),lim1,lim1).*dCphiM + repmat(Kb(ii,jj,:),lim1,lim1).*(dbound_greenM - CgradbpM);
%                     dLbeta(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Km(ii,jj,:),lim1,lim1).*dbeta;
%                 end
%             end
%             if N==1
%                 TX(lim1*ii-lim2:lim1*ii,s,:)=sum(dLX(lim1*ii-lim2:lim1*ii,:,:).*PHIX_s,2);
%                 if(SimType)
%                     TM(lim1*ii-lim2:lim1*ii,s,:)=sum(dLM(lim1*ii-lim2:lim1*ii,:,:).*PHIM_s,2);
%                     Tbeta(lim1*ii-lim2:lim1*ii,s,:)=sum(dLbeta(lim1*ii-lim2:lim1*ii,:,:).*PHIX_s,2);
%                 end
%             else if N==3
%                     TX(lim1*ii-lim2:lim1*ii,s,:)=sum(dLX(lim1*ii-lim2:lim1*ii,:,:).*[PHIX_s;PHIX_s],2);
%                     if(SimType)
%                         TM(lim1*ii-lim2:lim1*ii,s,:)=sum(dLM(lim1*ii-lim2:lim1*ii,:,:).*[PHIM_s;PHIM_s],2);
%                         Tbeta(lim1*ii-lim2:lim1*ii,s,:)=sum(dLbeta(lim1*ii-lim2:lim1*ii,:,:).*[PHIX_s;PHIX_s],2);
%                     end
%                 end
%             end
%         end
%
%         if (SimType)
%             delphi=-sum(PSIXMT.*repmat(TX(:,s,:),1,NofDet),1)+sum(PSIMMT.*repmat(Tbeta(:,s,:),1,NofDet),1)-sum(PSIMMT.*repmat(TM(:,s,:),1,NofDet),1);
%         else
%             delphi=-sum(PSIXXT.*repmat(TX(:,s,:),1,NofDet),1);
%         end
%         %         JacobianF(1:d,:)=delphi(1,1:d,:);
%
%         JacobianJ(1:NofDet,:)=delphi(1,1:NofDet,:);
%         JacJ=[JacJ;JacobianJ];
%         %         JacF=[JacF;JacobianF];
%     end
%
%     varargout{1}=JacJ;
%     % Fd=[Fd;FDet];
%     clear JacobianJ JacJ delphi CgradbpX
%     if (SimType)
%         clear CgradbpM dbeta
%     end
%     if mode=='H'
% %         tic
%         [d2CgradphiX,d2CgradbphiX,~]=d2SP3FormulationN(muax,musx,g_x,e_n,detE,Elem,nfibre,AOD); %Evaluate analytical second derivatives of the coeffecient matrices.
%          % Note that dCphi = zero.
%
%         if N>1
%             dICgradbphiX=get_dIC(CgradbphiX,muax,musx,g_x,e_n); % Use this to evaluate d/dx((Cgradb)^(-1)).
%         else
%             dICgradbphiX(:,:,e_n(:,1))=(1./CgradbphiX(:,:,e_n(:,1)))./(musx(:,:,e_n(:,1))*(1-g_x) + muax(:,:,e_n(:,1)));
%         end
%
%         Bound2_pertX=Cmat_zeros; % Boundary term in Term 1 of the SO derivative
%         Bound_dcba_X=Cmat_zeros; % Boundary term in SourceLambdaX - comes from perturbation eqn of the adjoint
%
%         for ii=1:size(e_n,1)
%             Bound2_pertX(:,:,e_n(ii,1))=+(d2CgradphiX(:,:,e_n(ii,1))/CgradbphiX(:,:,e_n(ii,1)))*CbphiX(:,:,e_n(ii,1)) ...
%                 - (((dCgradphiX(:,:,e_n(ii,1))/CgradbphiX(:,:,e_n(ii,1)))*dCgradbphiX(:,:,e_n(ii,1)))/CgradbphiX(:,:,e_n(ii,1)))*CbphiX(:,:,e_n(ii,1))...
%                 - ((CgradphiX(:,:,e_n(ii,1))*dICgradbphiX(:,:,e_n(ii,1))*dCgradbphiX(:,:,e_n(ii,1)))/CgradbphiX(:,:,e_n(ii,1)))*CbphiX(:,:,e_n(ii,1)) ...
%                 - (((CgradphiX(:,:,e_n(ii,1))/CgradbphiX(:,:,e_n(ii,1)))*d2CgradbphiX(:,:,e_n(ii,1)))/CgradbphiX(:,:,e_n(ii,1)))*CbphiX(:,:,e_n(ii,1)) ...
%                 - (CgradphiX(:,:,e_n(ii,1))/CgradbphiX(:,:,e_n(ii,1)))*dCgradbphiX(:,:,e_n(ii,1))*dICgradbphiX(:,:,e_n(ii,1))*CbphiX(:,:,e_n(ii,1)) ...
%                 + dCgradphiX(:,:,e_n(ii,1))*dICgradbphiX(:,:,e_n(ii,1))*CbphiX(:,:,e_n(ii,1));
%
%             Bound_dcba_X(:,:,e_n(ii,1))=(dCgradphiX(:,:,e_n(ii,1))/CgradbphiX(:,:,e_n(ii,1)))*CbphiX(:,:,e_n(ii,1)) ...
%                 + CgradphiX(:,:,e_n(ii,1))*dICgradbphiX(:,:,e_n(ii,1))*CbphiX(:,:,e_n(ii,1));
%         end
%
%         dLXpsi2=LZ;
%         for ii=1:D
%             for jj=1:D
%                 % for Lambda
%                 dLXpsi2(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)= repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiX + repmat(Km(ii,jj,:),lim1,lim1).*dCphiX + repmat(Kb(ii,jj,:),lim1,lim1).*(Bound_dcba_X);
%             end
%         end
%
%          %% Clear unnecessary variables. To save space
%             clear CgradphiX CphiX CgradbphiX CbphiX dCgradphiX dCgradbphiX dCphiX dICgradbphiX Bound_dcba_X dbound_greenX
%          %%
%
%         if (SimType)
%             [d2CgradphiM,d2CgradbphiM,~]=d2SP3FormulationN(muam,musm,g_m,e_n,detE,Elem,nfibre,AOD,muaf_mult^2);
%
%             if N>1
%                 dICgradbphiM=get_dIC(CgradbphiM,muam,musm,g_m,e_n,muaf_mult);
%             else
%                 dICgradbphiM(:,:,e_n(:,1))=(1./CgradbphiM(:,:,e_n(:,1)))./(musm(:,:,e_n(:,1))*(1-g_m) + muam(:,:,e_n(:,1)));
%             end
%             PSIXM=permute(conj(PSIXMT),[2,1,3]);
%             PSIMM=permute(conj(PSIMMT),[2,1,3]);
%
%             Bound2_pertM=Cmat_zeros;
%             Bound_dcba_M=Cmat_zeros;
%
%             for ii=1:size(e_n,1)
%                 Bound2_pertM(:,:,e_n(ii,1))=+(d2CgradphiM(:,:,e_n(ii,1))/CgradbphiM(:,:,e_n(ii,1)))*CbphiM(:,:,e_n(ii,1)) ...
%                     + dCgradphiM(:,:,e_n(ii,1))*dICgradbphiM(:,:,e_n(ii,1))*CbphiM(:,:,e_n(ii,1))...
%                     - (((dCgradphiM(:,:,e_n(ii,1))/CgradbphiM(:,:,e_n(ii,1)))*dCgradbphiM(:,:,e_n(ii,1)))/CgradbphiM(:,:,e_n(ii,1)))*CbphiM(:,:,e_n(ii,1))...
%                     - ((CgradphiM(:,:,e_n(ii,1))*dICgradbphiM(:,:,e_n(ii,1))*dCgradbphiM(:,:,e_n(ii,1)))/CgradbphiM(:,:,e_n(ii,1)))*CbphiM(:,:,e_n(ii,1)) ...
%                     - (((CgradphiM(:,:,e_n(ii,1))/CgradbphiM(:,:,e_n(ii,1)))*d2CgradbphiM(:,:,e_n(ii,1)))/CgradbphiM(:,:,e_n(ii,1)))*CbphiM(:,:,e_n(ii,1)) ...
%                     - (CgradphiM(:,:,e_n(ii,1))/CgradbphiM(:,:,e_n(ii,1)))*dCgradbphiM(:,:,e_n(ii,1))*dICgradbphiM(:,:,e_n(ii,1))*CbphiM(:,:,e_n(ii,1));
%
%                 Bound_dcba_M(:,:,e_n(ii,1))=(dCgradphiM(:,:,e_n(ii,1))/CgradbphiM(:,:,e_n(ii,1)))*CbphiM(:,:,e_n(ii,1)) ...
%                     + CgradphiM(:,:,e_n(ii,1))*dICgradbphiM(:,:,e_n(ii,1))*CbphiM(:,:,e_n(ii,1));
%
%             end
%
%             dLMpsi2=LZ;
%             for ii=1:D
%                 for jj=1:D
%                     % for Lambda
%                     dLMpsi2(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)= repmat(Ks(ii,jj,:),lim1,lim1).*dCgradphiM + repmat(Km(ii,jj,:),lim1,lim1).*dCphiM + repmat(Kb(ii,jj,:),lim1,lim1).*(Bound_dcba_M);
%                 end
%             end
%
%             TXT = zeros(((N+1)/2)*D,NofDet,Elem); TXlambda = TXT; TMT=TXT; TMlambda=TXT; Tbetalambda=TXT; % Initialisations
%             tempX=permute(conj(dLX),[2,1,3]); tempM=permute(conj(dLM),[2,1,3]);
%             temp2X=permute(conj(dLXpsi2),[2,1,3]);temp2M=permute(conj(dLMpsi2),[2,1,3]); temp3=permute(conj(dLbeta),[2,1,3]);
%             clear dLM dLX dLMpsi2 dLXpsi2 dLbeta LZ
%             for ii=1:D
%                 if N==1
%                     for dd=1:NofDet
%                         TXT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(tempX(lim1*ii-lim2:lim1*ii,:,:).*PSIXM(dd,:,:),2); % Source for thetaXM
%                         TXlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2X(lim1*ii-lim2:lim1*ii,:,:).*[PSIXM(dd,:,:)],2);
%                         TMT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(tempM(lim1*ii-lim2:lim1*ii,:,:).*PSIMM(dd,:,:),2); % Source for thetaMM
%                         TMlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2M(lim1*ii-lim2:lim1*ii,:,:).*[PSIMM(dd,:,:)],2);
%                         Tbetalambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp3(lim1*ii-lim2:lim1*ii,:,:).*[PSIMM(dd,:,:)],2); % beta term in source for theta XM
%                     end
%                 end
%
%                 if N==3
%                     for dd=1:NofDet
%                         %                     TXpsielem(lim1*ii-lim2:lim1*ii,dd,:)=-sum(dLXpsi.*[PSIXXTT(dd,:,:);PSIXXTT(dd,:,:)],2);
%                         TXT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(tempX(lim1*ii-lim2:lim1*ii,:,:).*[PSIXM(dd,:,:);PSIXM(dd,:,:)],2); % Source for thetaXM
%                         TXlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2X(lim1*ii-lim2:lim1*ii,:,:).*[PSIXM(dd,:,:);PSIXM(dd,:,:)],2);
%                         TMT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(tempM(lim1*ii-lim2:lim1*ii,:,:).*[PSIMM(dd,:,:);PSIMM(dd,:,:)],2);% Source for thetaMM
%                         TMlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2M(lim1*ii-lim2:lim1*ii,:,:).*[PSIMM(dd,:,:);PSIMM(dd,:,:)],2);
%                         Tbetalambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp3(lim1*ii-lim2:lim1*ii,:,:).*[PSIMM(dd,:,:);PSIMM(dd,:,:)],2); % beta term in source for thetaXM
%                     end
%                 end
%
%             end
%             SourceThetaXM=zeros(((N+1)/2)*Nodes,Elem,NofDet);
%             SourceLambdaX=zeros(((N+1)/2)*Nodes,Elem,NofSources);
%             SourceThetaMM=SourceThetaXM;
%             SourceLambdaM=SourceLambdaX;
%
%             %% Clear unnecessary variables. To save space
%             clear CgradphiM CphiM CgradbphiM CbphiM dCgradphiM dCgradbphiM dCphiM dICgradbphiM Bound_dcba_M dbound_greenM
%             clear tempX temp2X tempM temp2M temp3
%             %%
%
%         else
%             PSIXX=permute(conj(PSIXXT),[2,1,3]);
%             TXT = zeros(((N+1)/2)*D,NofDet,Elem); TXlambda = TXT;
%             temp=permute(conj(dLX),[2,1,3]);
%             temp2=permute(conj(dLXpsi2),[2,1,3]);
%             clear dLXpsi2 dLX
%             for ii=1:D
%                 if N==1
%                     for dd=1:NofDet
%                         TXT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(temp(lim1*ii-lim2:lim1*ii,:,:).*PSIXX(dd,:,:),2); % Source for thetaXX
%                         TXlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2(lim1*ii-lim2:lim1*ii,:,:).*[PSIXX(dd,:,:)],2);
%                     end
%                 end
%
%                 if N==3
%                     for dd=1:NofDet
%                         %                     TXpsielem(lim1*ii-lim2:lim1*ii,dd,:)=-sum(dLXpsi.*[PSIXXTT(dd,:,:);PSIXXTT(dd,:,:)],2);
%                         TXT(lim1*ii-lim2:lim1*ii,dd,:)=-sum(temp(lim1*ii-lim2:lim1*ii,:,:).*[PSIXX(dd,:,:);PSIXX(dd,:,:)],2); % Source for theta XX
%                         TXlambda(lim1*ii-lim2:lim1*ii,dd,:)=sum(temp2(lim1*ii-lim2:lim1*ii,:,:).*[PSIXX(dd,:,:);PSIXX(dd,:,:)],2);
%                     end
%                 end
%
%             end
%             clear temp temp2
%             SourceThetaXX=zeros(((N+1)/2)*Nodes,Elem,NofDet);
%             SourceLambdaX=zeros(((N+1)/2)*Nodes,Elem,NofSources);
%         end
%
%         for e=1:Elem
%             if ~(SimType)
%                 SourceThetaXX(DT(e,:),e,:)=permute(TXT(1:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%                 SourceLambdaX(DT(e,:),e,:)=-permute(TX(1:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%
%
%                 if N>=3
%                     SourceThetaXX(DT(e,:) + Nodes,e,:)=permute(TXT(2:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%                     SourceLambdaX(DT(e,:) + Nodes,e,:)=-permute(TX(2:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%                 end
%             else
%                 SourceThetaXM(DT(e,:),e,:)=permute((TXT(1:(N+1)/2:((N+1)/2)*D,:,e)+Tbetalambda(1:(N+1)/2:((N+1)/2)*D,:,e)),[1,3,2]);
%                 SourceLambdaX(DT(e,:),e,:)=-permute(TX(1:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%
%                 SourceThetaMM(DT(e,:),e,:)=permute(TMT(1:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%                 SourceLambdaM(DT(e,:),e,:)=-permute(TM(1:(N+1)/2:((N+1)/2)*D,:,e)-Tbeta(1:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%
%                 if N>=3
%                     SourceThetaXM(DT(e,:) + Nodes,e,:)=permute((TXT(2:(N+1)/2:((N+1)/2)*D,:,e)+Tbetalambda(2:(N+1)/2:((N+1)/2)*D,:,e)),[1,3,2]);
%                     SourceLambdaX(DT(e,:) + Nodes,e,:)=-permute(TX(2:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%
%                     SourceThetaMM(DT(e,:) + Nodes,e,:)=permute(TMT(2:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%                     SourceLambdaM(DT(e,:) + Nodes,e,:)=-permute(TM(2:(N+1)/2:((N+1)/2)*D,:,e)-Tbeta(2:(N+1)/2:((N+1)/2)*D,:,e),[1,3,2]);
%
%                 end
%
%             end
%         end
%         clear TXT TMT
%         LambdaX=zeros((N+1)/2*Nodes,Elem,NofSources);
%         parfor i=1:NofSources
%             LambdaX(:,:,i)=AX\(SourceLambdaX(:,:,i));
%         end
%
%         clear SourceLambdaX
%         if (SimType)
%             LambdaM=LambdaX;
%
%             parfor i=1:NofSources
%                 LambdaM(:,:,i)=AM\(SourceLambdaM(:,:,i))+AM\((MBetaT')*LambdaX(:,:,i));
%             end
%             clear SourceLambdaM
%             ThetaMM=zeros((N+1)/2*Nodes,Elem,NofDet);
%             ThetaXM=zeros((N+1)/2*Nodes,Elem,NofDet);
%             parfor i=1:NofDet
%                 ThetaMM(:,:,i)=(AM')\SourceThetaMM(:,:,i);
%             end
%             parfor i=1:NofDet
%                 ThetaXM(:,:,i)=(AX')\(SourceThetaXM(:,:,i))+(AX')\((MBetaT)*ThetaMM(:,:,i));
%             end
%               lambda_term_x=conj(TXlambda-Tbetalambda);
%               theta_term_m=TM-Tbeta;
%               clear TXlambda Tbetalambda TM Tbeta
%             clear SourceThetaXM SourceThetaMM MBetaT AX AM
%         else
%             ThetaXX=zeros((N+1)/2*Nodes,Elem,NofDet);
%             parfor i=1:NofDet
%                 ThetaXX(:,:,i)=AX'\(SourceThetaXX(:,:,i));
%             end
%             clear SourceThetaXX AX
%         end
%         %%
%
%         del2phi2=zeros(sum(1:Elem),NofDet*NofSources);
%         del2temp=zeros(sum(1:Elem),NofDet);
%         LZZ=zeros(((N+1)/2)*D,(N+1)/2*D,Elem);
%         TZZ=zeros((N+1)/2*D,1,Elem);
%         node1=T(:,1)-1;node2=T(:,2)-1; node3 =T(:,3)-1;
%         % Generate L matrices
%         LX2=LZZ; LM2 = LZZ;
%         clear LZZ
%          for ii=1:D
%                 for jj=1:D
%                     LX2(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*d2CgradphiX + repmat(Kb(ii,jj,:),lim1,lim1).*Bound2_pertX;
%                     if (SimType)
%                         LM2(lim1*ii-lim2:lim1*ii,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*d2CgradphiM + repmat(Kb(ii,jj,:),lim1,lim1).*Bound2_pertM;
%                     end
%                 end
%          end
%         clear d2CgradphiX Bound2_pertX muax0 muax musx edge_l PhiX PsiXX PSIXX Iv Jv CMat_zeros A
%         if SimType
%             clear d2CgradphiM Bound2_pertM muam muaxf musm PhiM PsiXM PsiMM PSIMM PSIXM Mbeta
%         end
%         for s=1:NofSources
%             PHIX_s=PHIX(s,:,:);
%             TX2elem=TZZ;
%             if (SimType)
%                 TM2elem=TX2elem;
%                 PHIM_s=PHIM(s,:,:);
%             end
%             for ii=1:D
% %                 for jj=1:D
% %                     LX2(:,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*d2CgradphiX + repmat(Kb(ii,jj,:),lim1,lim1).*Bound2_pertX;
% %                     if (SimType)
% %                         LM2(:,lim1*jj-lim2:lim1*jj,:)=repmat(Ks(ii,jj,:),lim1,lim1).*d2CgradphiM + repmat(Kb(ii,jj,:),lim1,lim1).*Bound2_pertM;
% %                     end
% %                 end
%                 if N==1
%                     TX2elem(lim1*ii-lim2:lim1*ii,1,:)=sum(LX2(lim1*ii-lim2:lim1*ii,:,:).*PHIX_s,2);
%                     if (SimType)
%                         TM2elem(lim1*ii-lim2:lim1*ii,1,:)=sum(LM2(lim1*ii-lim2:lim1*ii,:,:).*PHIM_s,2);
%                     end
%                 else if N==3
%                         TX2elem(lim1*ii-lim2:lim1*ii,1,:)=sum(LX2(lim1*ii-lim2:lim1*ii,:,:).*[PHIX_s;PHIX_s],2);
%                         if (SimType)
%                             TM2elem(lim1*ii-lim2:lim1*ii,1,:)=sum(LM2(lim1*ii-lim2:lim1*ii,:,:).*[PHIM_s;PHIM_s],2);
%                         end
%                     end
%                 end
%             end
%
%             if ~(SimType)
%                 d2T1X_temp = -sum(PSIXXT.*repmat(TX2elem,[1,NofDet,1]),1);
%                 temp_lambda=LambdaX(:,:,s);
%                 TX_temp=TX(:,s,:);
%
%                 parfor dd=1:NofDet
%                 del2temp(:,dd)=assemble_hess_vector(temp_lambda,conj(TXlambda(:,dd,:)),conj(ThetaXX(:,:,dd)),TX_temp,d2T1X_temp(:,dd,:),node1,node2,node3);
%                 end
%                 del2phi2(:,(s-1)*NofDet+1:s*NofDet) = del2temp;
% %                 clear del2temp;
% %                 TXelem_rep = repmat(TX(:,s,:),[1,d,1]);
%             else
%                 d2T1X_temp = -sum(PSIXMT.*repmat(TX2elem,[1,NofDet,1]),1);
%                 d2T1M_temp = -sum(PSIMMT.*repmat(TM2elem,[1,NofDet,1]),1);
%                 lambda_xs=LambdaX(:,:,s);
%                 lambda_ms=LambdaM(:,:,s);
%                 theta_term_x=TX(:,s,:);
%                 theta_term_ms=theta_term_m(:,s,:);
%
%                 parfor dd=1:NofDet
%                 del2temp(:,dd)=assemble_hess_vector(lambda_xs,lambda_term_x(:,dd,:),conj(ThetaXM(:,:,dd)),theta_term_x,d2T1X_temp(:,dd,:),node1,node2,node3);
%                 end
%                 del2phi2(:,(s-1)*NofDet+1:s*NofDet)=del2temp;
%                 parfor dd=1:NofDet
%                 del2temp(:,dd)=assemble_hess_vector(lambda_ms,conj(TMlambda(:,dd,:)),conj(ThetaMM(:,:,dd)),theta_term_ms,d2T1M_temp(:,dd,:),node1,node2,node3);
%                 end
%                 del2phi2(:,(s-1)*NofDet+1:s*NofDet)=del2phi2(:,(s-1)*NofDet+1:s*NofDet)+del2temp;
%
%             end
%         end
%             varargout{2}=del2phi2;
%             clear del2phi2 del2temp
%
%     end
%
% end
%
%
