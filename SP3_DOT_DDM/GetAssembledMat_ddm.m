function Tsparse=GetAssembledMat_ddm(sub_nodes,Dimension,N,Iv,Jv,Nodes,varargin)

% This function generates assembled matrices for SPN approximation 

%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub_nodes: Global node numbers of nodes lying in each subdomain
% Dimesnion - Choose 2 for 2D or 3 for 3D
% N  - order of the SPN approximation (3 for SP3 case)
% Iv,Jv: assembled matrices
%Nodes - number of nodes
% varargin: Variable input arguments contain elemental coefficient matrices (CgradphiX, CphiX and B1)
%           stiffness(Ks), mass(Km) and original boundary(Kb) assembled matrices

%C1,C2,C3 are elemental coefficient matrices in nodal basis
%K1,K2,K3 assembly matrices

%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%
% Tsparse: Assembled Sparse matrix for SPN formulation

% Ref 1 for better understanding of these assembled matrices for SP3 case
D=Dimension+1;
if nargin>=7
    C1=varargin{1};
    K1=varargin{2};
    T1=[];
    for ii=1:(N+1)/2
        T11=[];
        for jj=1:(N+1)/2
             T=repmat(C1(ii,jj,:),D,D).*K1;
             M=sparse(Iv(:),Jv(:),T(:),Nodes,Nodes);
             M=M(sub_nodes,sub_nodes);
             T11=[T11,M];
        end
        T1=[T1;T11];
    end
    Tsparse=T1;
end
if nargin>=9
    C2=varargin{3};
    K2=varargin{4};
    T2=[];
    for ii=1:(N+1)/2
        T21=[];
        for jj=1:(N+1)/2
             T=repmat(C2(ii,jj,:),D,D).*K2;
             M=sparse(Iv(:),Jv(:),T(:),Nodes,Nodes);
             M=M(sub_nodes,sub_nodes);
             T21=[T21,M];
        end
        T2=[T2;T21];
    end
    Tsparse=T1+T2;
end
if nargin>=11
    C3=varargin{5};
    K3=varargin{6};
    T3=[];
    for ii=1:(N+1)/2
        T31=[];
        for jj=1:(N+1)/2
             T=repmat(C3(ii,jj,:),D,D).*K3;
             M=sparse(Iv(:),Jv(:),T(:),Nodes,Nodes);
             M=M(sub_nodes,sub_nodes);
             T31=[T31,M];
        end
        T3=[T3;T31];
    end
    Tsparse=T1+T2+T3; % Tsparse = Ks*CgradphiX + Km*CphiX + Kb*B1 
end
