function dCgradbphi=Define_dBoundMat(N,n0,ni,mesh,dmu,mua1,varargin)

% Determine which edge numbers correspond to which boundary. 
dCgradbphi = zeros((N+1)/2,(N+1)/2,size(mesh.tri,1));

P = mesh.nodes;
if size(P,2)==2
    % 2-D case
    % Top
    top = ((P(mesh.edges(:,1),2)==min(mesh.nodes(:,2)))&(P(mesh.edges(:,2),2)==min(mesh.nodes(:,2))));
    % BElem = length(top);
    belem_top= mesh.edgeElem(top,1);
    right = (P(mesh.edges(:,1),1)==max(mesh.nodes(:,1))&(P(mesh.edges(:,2),1)==max(mesh.nodes(:,1))));
    belem_right = mesh.edgeElem(right,1);
    bottom = (P(mesh.edges(:,1),2)==max(mesh.nodes(:,2))&(P(mesh.edges(:,2),2)==max(mesh.nodes(:,2))));
    belem_bot = mesh.edgeElem(bottom,1);
    left = (P(mesh.edges(:,1),1)==min(mesh.nodes(:,1))&(P(mesh.edges(:,2),1)==min(mesh.nodes(:,1))));
    belem_left = mesh.edgeElem(left,1);
    switch(N)
        case 1
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1);
            
        
        case 3
            mua3 = varargin{1};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3);
            
        case 5
            mua3 = varargin{1};
            mua5 = varargin{2};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3,mua5);
            
        case 7
            mua3 = varargin{1};
            mua5 = varargin{2};
            mua7 = varargin{3};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3,mua5,mua7);
            
    end
else
    % 3D case
    % Top
    xmin = min(mesh.nodes(:,1)); ymin = min(mesh.nodes(:,2)); zmin = min(mesh.nodes(:,3));
    xmax = max(mesh.nodes(:,1)); ymax = max(mesh.nodes(:,2)); zmax = max(mesh.nodes(:,3));
    
    top = ((P(mesh.edges(:,1),3)==zmin)&(P(mesh.edges(:,2),3)==zmin)&(P(mesh.edges(:,3),3)==zmin));
    belem_top= mesh.edgeElem(top,1);
    right = ((P(mesh.edges(:,1),1)==xmax)&(P(mesh.edges(:,2),1)==xmax)&(P(mesh.edges(:,3),1)==xmax));
    belem_right= mesh.edgeElem(right,1);
    bottom = ((P(mesh.edges(:,1),3)==zmax)&(P(mesh.edges(:,2),3)==zmax)&(P(mesh.edges(:,3),3)==zmax));
    belem_bot= mesh.edgeElem(bottom,1);
    left = ((P(mesh.edges(:,1),1)==xmin)&(P(mesh.edges(:,2),1)==xmin)&(P(mesh.edges(:,3),1)==xmin));
    belem_left= mesh.edgeElem(left,1);
    front = ((P(mesh.edges(:,1),2)==ymin)&(P(mesh.edges(:,2),2)==ymin)&(P(mesh.edges(:,3),2)==ymin));
    belem_front= mesh.edgeElem(front,1);
    back= ((P(mesh.edges(:,1),2)==ymax)&(P(mesh.edges(:,2),2)==ymax)&(P(mesh.edges(:,3),2)==ymax));
    belem_back= mesh.edgeElem(back,1);
    
    switch(N)
        case 1
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1);
            [dCgradbphi(:,:,belem_front)] = fillmat(dmu,N,n0(5),ni,belem_front,mua1);
            [dCgradbphi(:,:,belem_back)] = fillmat(dmu,N,n0(6),ni,belem_back,mua1);
            
            
        case 3
            mua3 = varargin{1};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3);
            [dCgradbphi(:,:,belem_front)] = fillmat(dmu,N,n0(5),ni,belem_front,mua1,mua3);
            [dCgradbphi(:,:,belem_back)] = fillmat(dmu,N,n0(6),ni,belem_back,mua1,mua3);
            
        case 5
            mua3 = varargin{1};
            mua5 = varargin{2};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_front)] = fillmat(dmu,N,n0(5),ni,belem_front,mua1,mua3,mua5);
            [dCgradbphi(:,:,belem_back)] = fillmat(dmu,N,n0(6),ni,belem_back,mua1,mua3,mua5);
            
        case 7
            mua3 = varargin{1};
            mua5 = varargin{2};
            mua7 = varargin{3};
            [dCgradbphi(:,:,belem_top)] = fillmat(dmu,N,n0(1),ni,belem_top,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_right)] = fillmat(dmu,N,n0(2),ni,belem_right,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_bot)] = fillmat(dmu,N,n0(3),ni,belem_bot,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_left)] = fillmat(dmu,N,n0(4),ni,belem_left,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_front)] = fillmat(dmu,N,n0(5),ni,belem_front,mua1,mua3,mua5,mua7);
            [dCgradbphi(:,:,belem_back)] = fillmat(dmu,N,n0(6),ni,belem_back,mua1,mua3,mua5,mua7);
    end
end
end
function [dCgrad] = fillmat(dmu, N,n_out,n_in,belem, mua1,varargin)

nb = length(belem);
switch N
    case 1
        [~,B]=CalcReflCoeff(N,n_out,n_in);
        dCgrad=[-(1+B(1))./(3*mua1(belem).^2)].*dmu(belem);

        
    case 3
        mua3 = varargin{1};
        [~,B,~,D]=CalcReflCoeff(N,n_out,n_in);
        
        dCgrad=[-(1+B(1))*dmu(belem)./(3*mua1(belem).^2), D(1)*dmu(belem)./mua3(belem).^2
            D(2)*dmu(belem)./mua1(belem).^2, -(1+B(2))*dmu(belem)./(7*mua3(belem).^2)];
        
      
    case 5
        mua3 = varargin{1};
        mua5 = varargin{2};
        [~,B,~,D,~,F]=CalcReflCoeff(N,n_out,n_in);
        dCgrad=[-(1+B(1))*dmu(belem)./(3*mua1(belem).^2), D(1)*dmu(belem)./mua3(belem).^2,F(1)*dmu(belem)./mua5(belem).^2
            D(2)*dmu(belem)./mua1(belem).^2, -(1+B(2))*dmu(belem)./(7*mua3(belem).^2),F(2)./mua5(belem).^2
            D(3)*dmu(belem)./mua1(belem).^2,F(3)*dmu(belem)./mua3(belem).^2,-(1+B(3))*dmu(belem)./(11*mua5(belem).^2)];
        
     
    case 7
        mua3 = varargin{1};
        mua5 = varargin{2};
        mua7 = varargin{3};
        [~,B,~,D,~,F,~,H]=CalcReflCoeff(N,n_out,n_in); %Obtain boundary coefficients for SPN
        
        
        dCgrad=[-(1+B(1))*dmu(belem)./(3*mua1(belem).^2), D(1)*dmu(belem)./mua3(belem).^2,F(1)*dmu(belem)./mua5(belem).^2,H(1)*dmu(belem)./mua7(belem).^2
            D(2)*dmu(belem)./mua1(belem).^2, -(1+B(2))*dmu(belem)./(7*mua3(belem).^2),F(2)*dmu(belem)./mua5(belem).^2,H(2)*dmu(belem)./mua7(belem).^2
            D(3)*dmu(belem)./mua1(belem).^2,F(3)*dmu(belem)./mua3(belem).^2,-(1+B(3))*dmu(belem)./(11*mua5(belem).^2),H(3)*dmu(belem)./mua7(belem).^2
            D(4)*dmu(belem)./mua1(belem).^2,F(4)*dmu(belem)./mua3(belem).^2,H(4)*dmu(belem)./mua5(belem).^2,-(1+B(4))*dmu(belem)./(15*mua7(belem).^2)];
        
     
end
end