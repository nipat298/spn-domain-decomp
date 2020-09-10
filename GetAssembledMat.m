function Tsparse=GetAssembledMat(Dimension,N,Iv,Jv,Nodes,varargin)
%% generates assembled matrices 
%%%%%%%%%%%%%%%INPUT ARGUMENTS
%fname2- mat file cotnaining the assembly matrices
%Nodes - number of nodes
%C1,C2,C3 are matrices in elemental basis
%K1,K2,K3 assembly matrices
%%%%%%%%%%%%%OUTPUT ARGUMENTS
%Tsparse Sparse matrix.
%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
D=Dimension+1;
if nargin>=7
    C1=varargin{1};
    K1=varargin{2};
    T1=[];
    for ii=1:(N+1)/2
        T11=[];
        for jj=1:(N+1)/2
             T=repmat(C1(ii,jj,:),D,D).*K1;
             M=sparse(Iv,Jv,T(:),Nodes,Nodes);
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
             M=sparse(Iv,Jv,T(:),Nodes,Nodes);
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
             M=sparse(Iv,Jv,T(:),Nodes,Nodes);
             T31=[T31,M];
        end
        T3=[T3;T31];
    end
    Tsparse=T1+T2+T3;
end
