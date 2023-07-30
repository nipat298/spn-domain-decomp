function [Jd,varargout]=FwdSolver_DDM(Dimension,N,SimType,mesh,partmat,Src,Det,lambda,beta,AssemblyFile,options)
% This function implements the domain decomposition integrated forward solver for SP3 approximated DOT case.

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

  %%%%%% This case is specific to elastic scattering case SP3-DOT. The part corresponding to fluorescence 
  %%%%%% has been removed from here for sake of clarity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
else
    % Elastic scattering case
    
    % Obtaining the optical coefficient value for all elements
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
    

    %%%% The following commented lines elaborate 'assembly on fly' mechanism which was used, in some cases, 
    %%%% by Nishigandha while developing codes for SPN approximation
    %%%% This approach has not been used for domain decomposition codes. Can be explored in cases of severe 
    %%%% memory constraints where storing the large sized assembled matrix could be cumbersome.

    %%%% if (isa(AssemblyFile,'function_handle'))
    %%%%         disp('Obtain FEM matrices...')
    %%%%         [Ks,Km,Kb,Iv,Jv,~] = AssemblyFile(mesh);
    %%%%     else
    %%%%         load(AssemblyFile,'Ks','Km','Kb','Iv','Jv');
    %%%%     end
    
   
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
        
        disp('Generate Sparse Assembly matrix... ')
        sub_nodes=partmat(ii).nodes; % Global node numbers of nodes lying in each subdomain
        partmat(ii).Gsub = GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,CgradphiX(:,:,elem_sub),Ks,CphiX(:,:,elem_sub),Km,B1,Kb);
        % Gsub stores the assembled Tsparse matrix for each subdomain

        % This function obtains assembly matrix for interface nodes for each subdomain.
        % This matrix is added both to the LHS and RHS of the main equation.That part is coovered in the following two sections.
        ABC=interface_mm(partmat(ii),mesh.nodes,NoDomains,beta);          
         for j=1:NoDomains
             temp_int=ABC.RHS_extra{j};
              % Obtaining SP3 format compatible matrix for the assembled interface matrix for RHS
             temp_mat=[temp_int, zeros(size(temp_int));zeros(size(temp_int)), temp_int];
             partmat(ii).RHS_extra{j}=temp_mat; % Adding the assembled interface matrix to RHS 
         end
         
         % Adding the assembled interface matrix to LHS
          node_len=length(sub_nodes);
          LHS_extra=sparse(node_len,node_len);
          % Obtaining SP3 format compatible matrix for the assembled interface matrix for LHS
          LHS_extra=[LHS_extra, zeros(size(LHS_extra));zeros(size(LHS_extra)), LHS_extra];
        for j=1:NoDomains
            if j==ii
                continue
            else
                LHS_extra=LHS_extra+partmat(ii).RHS_extra{j};  
            end
        end
        partmat(ii).LHS=(lambda)*LHS_extra; %Interface matrix is added to LHS with a multiplication of lambda
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
    % step of assembling matrices for each subdomains.
    
    % This involves implementing the fundamental part of domain decomposition algorithm 
    partmat=OptimalSchwarz(partmat,lambda,N);
    
    % The following section helps obtain the fluence values for SP3 case.
    % The SP3 solution vector is twice the size of the solution vector to incorporate results for both 1st and 2nd composite moment
    PhiX=zeros(2*Nodes,1);
    for ii=1:NoDomains
    PhiX(partmat(ii).nodes)=partmat(ii).soln(1:length(partmat(ii).nodes));
    
    n11(:)=partmat(ii).nodes(:)+ Nodes;
    n11=n11';
    PhiX(n11)=partmat(ii).soln((1+length(partmat(ii).nodes)):end);
    end

% The above section was the final part of DDM algorithm.
% We obtain the fluence values for domain decomposition integrated SP3 DOT forward solver.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
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
