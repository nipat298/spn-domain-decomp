function [Jd,varargout]=FwdSolver_DDM(Dimension,N,SimType,mesh,partmat,Src,Det,lambda,beta,AssemblyFile,options)
% This function implements the domain decomposition integrated forward solver for SP3 approximated FOT case.

%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%
% Dimesnion - Choose 2 for 2D or 3 for 3D
% N  - order of the SPN approximation (3 for SP3 case)
% mesh - structure type variable containing meshing information
% Src - structure type variable containing the source information
% Det - structure type variable containing the detector information
% lambda - transmission coefficient at the interface
% beta - penalty coefficient (extra parameter included by Deng-method 2)
% AssemblyFile - name of the mat file containing the mass,stiffness matrices and indexing vectors
% evalFluence - Flag set to '0' to suppress fluence output, else set to '1' to
% probType - 'fwd' for data generation, forward solve only,
%            'recon' for generating data in a reconstruction routine. saves some variables that are 
%             needed in the reconstruction routine as well as evaluation of adjoint sensitivity

%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%
% Jd - Exiting partial current if evalFluence ==1
% varargout{1} = Fluence
% varargout{2} = FwData (if probType = 'recon') - structure containing information required for sensitivity and reconstruction routine
% if evalFluence == 0 and probType = 'recon'
% varargout{1} = FwData
%%

mesh_spacing=mean([mesh.hx,mesh.hy,mesh.hz]);
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

    % Coefficient matrices
    [CgradphiX,CphiX,CgradbphiX,CbphiX,CJX]=SP3Formulation(N,Elem,mesh,muax,mesh.opt.musx,bn,mesh.opt.refrind_out,mesh.opt.refrind,Src,Det); %Obtain the SP3 matrices for excitation
    
    NoDomains=size(partmat,2); % Refers to number of subdomains
    
        for ii=1:NoDomains
        [Ks,Km,Kb,Iv,Jv]=assemble_ddm(mesh,partmat,ii); % Assembly of stiffness,mass and boundary matrices for each subdomain
        elem_sub=partmat(ii).ECM(:,5); % Global numbering of elements in each subdomain
        elem_sublen=length(partmat(ii).ECM(:,5));
        B=zeros((N+1)/2,(N+1)/2,elem_sublen);    B1=B;

        % The next two lines give the elements lying along the original boundary of a subdomain
        ogbound_elem=partmat(ii).belem(:,1);
        ogbound_elemGB=partmat(ii).ECM(ogbound_elem,5);

        % The following section calculates elemental coefficient matrix for the boundary matrix (original subdomain boundary). 
        % Ref 1 for better understanding of these equations
        for i=1:size(ogbound_elem,1)
            B(:,:,ogbound_elem(i,1))=CgradphiX(:,:,ogbound_elem(i,1))/CgradbphiX(:,:,ogbound_elemGB(i,1));
            B1(:,:,ogbound_elem(i,1))=B(:,:,ogbound_elem(i,1))*CbphiX(:,:,ogbound_elemGB(i,1));
        end
        
        disp('Generate Sparse Assembly matrix for excitation... ')
        sub_nodes=partmat(ii).nodes;  % Global node numbers of nodes lying in each subdomain
        partmat(ii).Gsub = GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,CgradphiX(:,:,elem_sub),Ks,CphiX(:,:,elem_sub),Km,B1,Kb);
        % Gsub stores the assembled Tsparse matrix for each subdomain


        % This function obtains assembly matrix for interface nodes for each subdomain.
        % This matrix is added both to the LHS and RHS of the main equation.That part is coovered in the following two sections.
        ABC=interface_mm(partmat(ii),mesh.nodes,NoDomains,beta);
         
         for j=1:NoDomains
             temp_int=ABC.RHS_extra{j};
             ss=size(temp_int,1);
             % Obtaining SP3 format compatible matrix for the assembled interface matrix for RHS
             temp_mat=[temp_int, sparse(ss,ss);sparse(ss,ss), temp_int];
             partmat(ii).RHS_extra{j}=temp_mat;  % Adding the assembled interface matrix to RHS 
         end
         
         % Adding the assembled interface matrix to LHS
          node_len=length(sub_nodes);
          LHS_extra=sparse(node_len,node_len);
          ss=size(LHS_extra,1);
          % Obtaining SP3 format compatible matrix for the assembled interface matrix for LHS
          LHS_extra=[LHS_extra, sparse(ss,ss);sparse(ss,ss), LHS_extra];
        for j=1:NoDomains
            if j==ii
                continue
            else
                LHS_extra=LHS_extra+partmat(ii).RHS_extra{j};
            end
        end
        partmat(ii).LHS=(lambda)*LHS_extra;%Interface matrix is added to LHS with a multiplication of lambda
        partmat(ii).Gsub=partmat(ii).Gsub+partmat(ii).LHS;
        % Gsub is the final matrix on LHS which is formed by adding all matrices on LHS        
         
        % The following section performs source calculation
        for s=1:NofSources
           elemindex=ismember(partmat(ii).ECM(:,5),Src(s).elem);
            if ~isempty(elemindex)
                sub_nodes=partmat(ii).nodes;
                SX_sub=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,Src(s));
                partmat(ii).source(:,s)=SX_sub([sub_nodes;sub_nodes+Nodes]);
            end
        end
        
        end % This ends the loop for number of subdomains and concludes the pre processing
    % step of assembling matrices (for excitation equations) for each subdomains.
    
    
    % This involves implementing the fundamental part of domain decomposition algorithm 
    % (for excitation equation)
    % OptimalSchwarz_ex: this subroutine computes solution for excitation equations through DDM  
    partmat=OptimalSchwarz_ex(partmat,lambda,N);
    
      % The following section helps obtain the fluence values for excitation equation of SP3 FOT case.
    % The SP3 solution vector is twice the size of the solution vector to incorporate results for both 1st and 2nd composite moment
    PhiX=zeros(2*Nodes,1);
    for ii=1:NoDomains
    PhiX(partmat(ii).nodes)=partmat(ii).soln(1:length(partmat(ii).nodes));
    
    n11(:)=partmat(ii).nodes(:)+ Nodes;
    n11=n11';
    PhiX(n11)=partmat(ii).soln((1+length(partmat(ii).nodes)):end);
   partmat(ii).Phi_emission=partmat(ii).soln;
   end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%This concludes the DDM incorporated solution of excitation equation for SP3 FOT case%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%% From here begins the solution for emission equation %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    for ii=1:NoDomains
        % Assembly of stiffness,mass and boundary matrices for emission for each subdomain
        [Ks,Km,Kb,Iv,Jv]=assemble_ddm(mesh,partmat,ii);

        elem_sub=partmat(ii).ECM(:,5);
        elem_sublen=length(partmat(ii).ECM(:,5));
        B=zeros((N+1)/2,(N+1)/2,elem_sublen);    B1=B;
        %Generates the boundary term in equation (57)
        ogbound_elem=partmat(ii).belem(:,1);
        ogbound_elemGB=partmat(ii).ECM(ogbound_elem,5);
        
        % These are the same steps as were followed for excitation equation
        for i=1:size(ogbound_elem,1)
            B(:,:,ogbound_elem(i,1))=CgradphiM(:,:,ogbound_elem(i,1))/CgradbphiM(:,:,ogbound_elemGB(i,1));
            B1(:,:,ogbound_elem(i,1))=B(:,:,ogbound_elem(i,1))*CbphiM(:,:,ogbound_elemGB(i,1));
        end
    
     disp('Generate Sparse Assembly matrix for emission... ')
        sub_nodes=partmat(ii).nodes;
        partmat(ii).Gsub = GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,CgradphiM(:,:,elem_sub),Ks,CphiM(:,:,elem_sub),Km,B1,Kb);
        partmat(ii).Gsub=partmat(ii).Gsub+partmat(ii).LHS;
        partmat(ii).SM=GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,Mbeta(:,:,elem_sub),Km)*partmat(ii).Phi_emission;
      % partmat(ii).SM obtains the RHS term for emission equation. 
      % The solution of the ecxitation equation serves as the source for emission equation
    end
    
    % The interface matrix that was calculated through "interface_mm" subroutine is 
    % already stored in the 'partmat' structure. Therefore, we are not calculating it again for emission equation. 
    
    % OptimalSchwarz_em: this subroutine computes solution for emission equations through DDM  
    partmat=OptimalSchwarz_em(partmat,lambda,N);
    
    PhiM=zeros(2*Nodes,1);
    for ii=1:NoDomains
    PhiM(partmat(ii).nodes)=partmat(ii).solnFOT(1:length(partmat(ii).nodes));
    
    n11(:)=partmat(ii).nodes(:)+ Nodes;
    n11=n11';
    PhiM(n11)=partmat(ii).solnFOT((1+length(partmat(ii).nodes)):end);
    end
    
% The above section was the concluding step of SP3 FOT DDM algorithm.
% We obtain the fluence values for emission equation at this step.
   
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
    
%     if (isa(AssemblyFile,'function_handle'))
%         disp('Obtain FEM matrices...')
%         [Ks,Km,Kb,Iv,Jv,~] = AssemblyFile(mesh);
%     else
%         load(AssemblyFile,'Ks','Km','Kb','Iv','Jv');
%     end
    
     NoDomains=size(partmat,2);
    
%     beta=1/(mesh_spacing^2);
%     alpha=1/(mesh_spacing);
 
    for ii=1:NoDomains
        [Ks,Km,Kb,Iv,Jv]=assemble_ddm(mesh,partmat,ii);
        elem_sub=partmat(ii).ECM(:,5);
        elem_sublen=length(partmat(ii).ECM(:,5));
        B=zeros((N+1)/2,(N+1)/2,elem_sublen);    B1=B;
        %Generates the boundary term in equation (57)
        ogbound_elem=partmat(ii).belem(:,1);
        ogbound_elemGB=partmat(ii).ECM(ogbound_elem,5);
        for i=1:size(ogbound_elem,1)
            B(:,:,ogbound_elem(i,1))=CgradphiX(:,:,ogbound_elem(i,1))/CgradbphiX(:,:,ogbound_elemGB(i,1));
            B1(:,:,ogbound_elem(i,1))=B(:,:,ogbound_elem(i,1))*CbphiX(:,:,ogbound_elemGB(i,1));
        end
        
        disp('Generate Sparse Assembly matrix... ')
        sub_nodes=partmat(ii).nodes;
        partmat(ii).Gsub = GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,CgradphiX(:,:,elem_sub),Ks,CphiX(:,:,elem_sub),Km,B1,Kb);
        
%         temp_int=interface_mm(partmat(ii),mesh.nodes,NoDomains);
% % % %          temp_int=interface_mm(partmat(ii),mesh.nodes,NoDomains);
% % % %          for j=1:NoDomains
% % % %             partmat(ii).RHS_extra(j)=[temp_int(j), zeros(size(temp_int(j)));zeros(size(temp_int(j))), temp_int(j)];
% % % %         end
         %partmat(ii).RHS_extra=[temp_int, zeros(size(temp_int));zeros(size(temp_int)), temp_int];
         
         ABC=interface_mm(partmat(ii),mesh.nodes,NoDomains,beta);
         
         for j=1:NoDomains
             temp_int=ABC.RHS_extra{j};
             temp_mat=[temp_int, zeros(size(temp_int));zeros(size(temp_int)), temp_int];
             partmat(ii).RHS_extra{j}=temp_mat;
%              partmat(ii).RHS_extra(j)=[temp_int, zeros(size(temp_int));zeros(size(temp_int)), temp_int];
         end
         
         
          node_len=length(sub_nodes);
          LHS_extra=sparse(node_len,node_len);
          LHS_extra=[LHS_extra, zeros(size(LHS_extra));zeros(size(LHS_extra)), LHS_extra];
        for j=1:NoDomains
            if j==ii
                continue
            else
                LHS_extra=LHS_extra+partmat(ii).RHS_extra{j};
            end
        end
        partmat(ii).LHS=(lambda)*LHS_extra;
        partmat(ii).Gsub=partmat(ii).Gsub+partmat(ii).LHS;
%         partmat(ii).Gsub=A_sub(sub_nodes,sub_nodes);
        
%         for jj=1:size(partmat,1)
%              ismem=ismember(partmat(i).art_bound,partmat(i).interface_g{j});
%             indx=(ismem(:,1)==1 & ismem(:,2)==1 & ismem(:,3)==1);
%             temp=partmat(i).art_bound(indx,:);
%             
%        art_int_sparse=GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,imat,Kb_art);
%         
         
        
        for s=1:NofSources
%             elemindex=find(partmat(ii).ECM(:,5)==Src(s).elem);
            elemindex=ismember(partmat(ii).ECM(:,5),Src(s).elem);
            if ~isempty(elemindex)
                sub_nodes=partmat(ii).nodes;
                SX_sub=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,Src(s));
                partmat(ii).source(:,s)=SX_sub([sub_nodes;sub_nodes+Nodes]);
            end
        end
        
    end
    
    
% load('partmat_OS')
    
    partmat=OptimalSchwarz(partmat,lambda,N);
    
    % 
    PhiX=zeros(2*Nodes,1);
    for ii=1:NoDomains
    PhiX(partmat(ii).nodes)=partmat(ii).soln(1:length(partmat(ii).nodes));
    
    n11(:)=partmat(ii).nodes(:)+ Nodes;
    n11=n11';
    PhiX(n11)=partmat(ii).soln((1+length(partmat(ii).nodes)):end);
    end

% PhiX(partmat(2).nodes)=partmat(2).soln(1:length(partmat(2).nodes)); %Because nodes from second subdomain are being assigned value after first, so interface nodes would automatically take up the fluence values from 2nd subdomain
% n22(:)=partmat(2).nodes(:)+ Nodes;
% n22=n22';
% PhiX(n22)=partmat(2).soln((1+length(partmat(2).nodes)):end);
% 
% PhiX(partmat(3).nodes)=partmat(3).soln(1:length(partmat(3).nodes));
% n33(:)=partmat(3).nodes(:)+ Nodes;
% n33=n33';
% PhiX(n33)=partmat(3).soln((1+length(partmat(3).nodes)):end);
% 
% PhiX(partmat(4).nodes)=partmat(4).soln(1:length(partmat(4).nodes));
% n44(:)=partmat(4).nodes(:)+ Nodes;
% n44=n44';
% PhiX(n44)=partmat(4).soln((1+length(partmat(4).nodes)):end);
%     
    
    %Solve the forward problem to compute PhiX and PhiM
    
    %AX=Spnfwd(Dimension,N,CgradphiX,CphiX,CgradbphiX,CbphiX,mesh,Iv,Jv,Ks,Km,Kb);
    
    %     SX=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,mesh.edgeElem,mesh.edgeType,Src);
%     for s=1:NofSources
%         SX(:,s)=GetSource(Dimension,N,CgradphiX,CgradbphiX,mesh.nodes,mesh.tri,Src(s));
%         %     PhiX(:,s) = gmres(AX,SX(:,s),gm_restart, gm_tol, gm_maxit);
%     end
%     
%     switch options.itersolver
%         case 'none'
%             PhiX = AX\SX;
%         case 'gmres'
%             for s=1:NofSources
%                 PhiX(:,s) = gmres(AX,SX(:,s),options.restart,options.tol, options.maxit);
%             end
%         case 'bicgstab'
%             for s=1:NofSources
%                 PhiX(:,s) = bicgstab(AX,SX(:,s),options.tol, options.maxit);
%             end
%         case 'pcg'
%             for s=1:NofSources
%                 PhiX(:,s) = pcg(AX,SX(:,s),options.tol, options.maxit);
%             end
%     end
%     
    
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
