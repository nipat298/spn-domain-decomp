function [Jd,varargout]=FwdSolver(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,options)
%%
% Modified on 4th Dcember 2020 to support extended sources
% Modified on 9th April 2019
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
    
    switch(options.pfuncflag)
        case 0
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx);
        case 1
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.gx, mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.gx,mesh.opt.alpha,mesh.opt.gbx);
      
    end
    [CgradphiX,CphiX,CgradbphiX,CbphiX,CJX]=SP3Formulation(N,Elem,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for excitation
    
     if (isa(AssemblyFile,'function_handle'))
        disp('Obtain FEM matrices...')
        [Ks,Km,Kb,Iv,Jv] = AssemblyFile(mesh);
    else
        load(AssemblyFile,'Ks','Km','Kb','Iv','Jv');
    end
   
    
    
    %Solve the forward problem to compute PhiX and PhiM
    AX=Spnfwd(Dimension,N,CgradphiX,CphiX,CgradbphiX,CbphiX,mesh,Iv,Jv,Ks,Km,Kb);
    
    for s=1:NofSources
           SX(:,s)=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,Src(s));        
    end
%     SX = ones(2*Nodes,1);
%     PhiX=AX\SX;
switch options.itersolver
    case 'none'
        PhiX = AX\SX; 
    case 'gmres'
        for s=1:NofSources
            PhiX(:,s) = gmres(AX,SX(:,s),options.restart,options.tol, options.maxit);
        end
    case 'bicgstab'
        for s=1:NofSources
            PhiX(:,s) = bicgstab(AX,SX(:,s),options.tol, options.maxit);
        end
    case 'pcg'
        for s=1:NofSources
            PhiX(:,s) = pcg(AX,SX(:,s),options.tol, options.maxit);
        end
end
    % Emission field
    muam=zeros(1,1,Elem); %absorption coefficient at emission wavelength
    muam(1,1,:)=mesh.opt.muaimult*mesh.opt.muaxi(1,1,:)+mesh.opt.muafmult*mesh.opt.muaxf(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,Elem);   
    
switch(options.pfuncflag)
    case 0
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm);
    case 1
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm,mesh.opt.gm,mesh.opt.gamma);
    case 2
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gm, mesh.opt.gm,mesh.opt.alpha, mesh.opt.gbm);
end
    [CgradphiM,CphiM,CgradbphiM,CbphiM,CJM]=SP3Formulation(N,Elem,mesh,muam,mesh.opt.musm,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for emission
    mesh.opt.beta=(mesh.opt.Qf/(1+1i*Src(1).mfreq*mesh.opt.tau)).*mesh.opt.muaxf(1,1,:); %(muax0(1,1,:) + muaxf(1,1,:)); %Note this is muaxf and not muax as in earlier versions
    switch(N)
        case 1
            Mbeta = mesh.opt.beta;
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
    
    AM=Spnfwd(Dimension,N,CgradphiM,CphiM,CgradbphiM,CbphiM,mesh,Iv,Jv,Ks,Km,Kb);
    % Obtain emission source
    SM = GetAssembledMat(Dimension,N,Iv,Jv,Nodes,Mbeta,Km)*PhiX;
    
    switch options.itersolver
        case 'none'
            PhiM = AM\SM;
        case 'gmres'
            for s=1:NofSources
                PhiM(:,s) = gmres(AM,SM(:,s),options.restart,options.tol, options.maxit);
            end
        case 'bicgstab'
            for s=1:NofSources
                PhiM(:,s) = bicgstab(AM,SM(:,s),options.tol, options.maxit);
            end
        case 'pcg'
            for s=1:NofSources
                PhiM(:,s) = pcg(AM,SM(:,s),options.tol, options.maxit);
            end
    end
    
    
    if (options.evalFluence)
        switch(N)
            case 1
                FX = PhiX(1:Nodes,:); 
                FM = PhiM(1:Nodes,:);
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
        JmatX=GetPartCurr(Dimension,N,CJX,mesh.nodes,T,Det);
        JX=JmatX*PhiX;
        JmatM=GetPartCurr(Dimension,N,CJM,mesh.nodes,T,Det);
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
      
     if strcmp(options.probType,'recon')
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
        if (options.evalFluence)
          varargout{2} = FwData;
        else
            varargout{1} = FwData;
        end
    end

    
    
else
    % Elastic scattering case
    
    muax=zeros(1,1,Elem);
    muax(1,1,:)=mesh.opt.muaxi(1,1,:)+(1i*Src(1).mfreq/mesh.opt.speed)*ones(1,1,Elem);  
    switch(options.pfuncflag)
        case 0
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx);
        case 1
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.gx,mesh.opt.gamma);
        case 2
            bn = get_pfunc_coeff(options.pfuncflag,deflag,mesh.opt.gx,mesh.opt.gx,mesh.opt.alpha, mesh.opt.hx);
    end
    % Coefficient matrices
    [CgradphiX,CphiX,CgradbphiX,CbphiX,CJX]=SP3Formulation(N,Elem,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for excitation
   
     if (isa(AssemblyFile,'function_handle'))
        disp('Obtain FEM matrices...')
        [Ks,Km,Kb,Iv,Jv,~] = AssemblyFile(mesh);
    else
        load(AssemblyFile,'Ks','Km','Kb','Iv','Jv');
    end
   
    %Solve the forward problem to compute PhiX and PhiM
    
    AX=Spnfwd(Dimension,N,CgradphiX,CphiX,CgradbphiX,CbphiX,mesh,Iv,Jv,Ks,Km,Kb);
    
%     SX=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,mesh.edgeType,Src);
    for s=1:NofSources
    SX(:,s)=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,Src(s));
%     PhiX(:,s) = gmres(AX,SX(:,s),gm_restart, gm_tol, gm_maxit); 
    end
  
    switch options.itersolver
        case 'none'
            PhiX = AX\SX;
        case 'gmres'
            for s=1:NofSources
                PhiX(:,s) = gmres(AX,SX(:,s),options.restart,options.tol, options.maxit);
            end
        case 'bicgstab'
            for s=1:NofSources
                PhiX(:,s) = bicgstab(AX,SX(:,s),options.tol, options.maxit);
            end
        case 'pcg'
            for s=1:NofSources
                PhiX(:,s) = pcg(AX,SX(:,s),options.tol, options.maxit);
            end
    end
    
    
%     PhiX=AX\SX;
%     PhiX = gmres(AX,SX,gm_restart, gm_tol, gm_maxit); 
    
%     JmatX=GetPartCurr(Dimension,N,CJX,mesh.nodes,T,mesh.edgeElem,mesh.edgeType,Det);
    JmatX=GetPartCurr(Dimension,N,CJX,mesh.nodes,T,Det);
    JX =JmatX*PhiX;
%     Jd = JX(:);    
     if (options.evalFluence)
        switch(N)
            case 1 
                FX = PhiX(1:Nodes,:);
            case 3
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:end,:);            
            case 5
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:end,:);             
            case 7
                FX = PhiX(1:Nodes,:) - (2/3)*PhiX(Nodes+1:2*Nodes,:) + (8/15)*PhiX(2*Nodes+1:3*Nodes) - (16/35)*PhiX(3*nodes+1:end,:);
        end 
     end
    Jd = cell(NofSources,1);
for s=1:NofSources
    Jd{s} = JX(:,s);
end
% Jd = JX; 
     if (options.evalFluence)
        varargout{1}=PhiX; %FX;
        
     end
     
     if strcmp(options.probType,'recon')
        FwData.cgradphix= CgradphiX;
        FwData.cphix= CphiX;
        FwData.cgradbx= CgradbphiX;
        FwData.cbx= CbphiX;
        FwData.cjx = CJX;
        FwData.sparsex = AX;
        FwData.phix = PhiX;
         if (options.evalFluence)
          varargout{2} = FwData;
        else
            varargout{1} = FwData;
        end
    end

end
fprintf('\nDone... \n');
