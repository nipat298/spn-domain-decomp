function [Cgradbphi,Cbphi]=DefineBoundMat(N,n0,ni,mesh,mua1,varargin)
%%
% Modified on 22 December to support variable N
% function to Determine which edge numbers correspond to which boundary.
% see SP3Formulation for definitions of input and output arguments
%%
Cgradbphi = zeros((N+1)/2,(N+1)/2,size(mesh.tri,1));
Cbphi = Cgradbphi;
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
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1);
        case 3
             mua3 = varargin{1};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3);
            
        case 5
             mua3 = varargin{1};
            mua5 = varargin{2};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3,mua5);
            
        case 7
             mua3 = varargin{1};
            mua5 = varargin{2};
            mua7 = varargin{3};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3,mua5,mua7);
            
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
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1);
            [Cgradbphi(:,:,belem_front), Cbphi(:,:,belem_front)] = fillmat(N,n0(5),ni,belem_front,mua1);
            [Cgradbphi(:,:,belem_back), Cbphi(:,:,belem_back)] = fillmat(N,n0(6),ni,belem_back,mua1);
            
        case 3
             mua3 = varargin{1};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3);
            [Cgradbphi(:,:,belem_front), Cbphi(:,:,belem_front)] = fillmat(N,n0(5),ni,belem_front,mua1,mua3);
            [Cgradbphi(:,:,belem_back), Cbphi(:,:,belem_back)] = fillmat(N,n0(6),ni,belem_back,mua1,mua3);
            
        case 5
            mua3 = varargin{1};
            mua5 = varargin{2};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_front), Cbphi(:,:,belem_front)] = fillmat(N,n0(5),ni,belem_front,mua1,mua3,mua5);
            [Cgradbphi(:,:,belem_back), Cbphi(:,:,belem_back)] = fillmat(N,n0(6),ni,belem_back,mua1,mua3,mua5);
            
        case 7
              mua3 = varargin{1};
            mua5 = varargin{2};
            mua7 = varargin{3};
            [Cgradbphi(:,:,belem_top), Cbphi(:,:,belem_top)] = fillmat(N,n0(1),ni,belem_top,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_right), Cbphi(:,:,belem_right)] = fillmat(N,n0(2),ni,belem_right,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_bot), Cbphi(:,:,belem_bot)] = fillmat(N,n0(3),ni,belem_bot,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_left), Cbphi(:,:,belem_left)] = fillmat(N,n0(4),ni,belem_left,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_front), Cbphi(:,:,belem_front)] = fillmat(N,n0(5),ni,belem_front,mua1,mua3,mua5,mua7);
            [Cgradbphi(:,:,belem_back), Cbphi(:,:,belem_back)] = fillmat(N,n0(6),ni,belem_back,mua1,mua3,mua5,mua7);
    end
end
end
function [Cgrad, C] = fillmat(N,n_out,n_in,belem, mua1,varargin)
nb = length(belem);
switch N
    case 1
        [A,B]=CalcReflCoeff(N,n_out,n_in);
        Cgrad=(1+B(1))./(3*mua1(belem));
        C = (0.5+A(1))*ones(1,1,nb); 
    case 3
             mua3 = varargin{1};
        [A,B,C,D]=CalcReflCoeff(N,n_out,n_in);
        
        Cgrad=[(1+B(1))./(3*mua1(belem)), -D(1)./mua3(belem)
            -D(2)./mua1(belem), (1+B(2))./(7*mua3(belem))];
        
        C=[(0.5+A(1))*ones(1,1,nb), -(1/8+C(1))*ones(1,1,nb)
            -(1/8+C(2))*ones(1,1,nb),(7/24+A(2))*ones(1,1,nb)];
    case 5
       mua3 = varargin{1};
            mua5 = varargin{2};
        [A,B,C,D,E,F]=CalcReflCoeff(N,n_out,n_in);
        Cgrad=[(1+B(1))./(3*mua1(belem)), -D(1)./mua3(belem),-F(1)./mua5(belem)
            -D(2)./mua1(belem), (1+B(2))./(7*mua3(belem)),-F(2)./mua5(belem)
            -D(3)./mua1(belem),-F(3)./mua3(belem),(1+B(3))./(11*mua5(belem))];
        
        C=[(0.5+A(1))*ones(1,1,nb), -(1/8+C(1))*ones(1,1,nb),-(-1/16+E(1))*ones(1,1,nb)
            -(1/8+C(2))*ones(1,1,nb),(7/24+A(2))*ones(1,1,nb),-(41/384+E(2))*ones(1,1,nb)
            -(-1/16+C(3))*ones(1,1,nb),-(41/384+E(3))*ones(1,1,nb),(407/1920+A(3))*ones(1,1,nb)];
    case 7
          mua3 = varargin{1};
            mua5 = varargin{2};
            mua7 = varargin{3};
        [A,B,C,D,E,F,G,H]=CalcReflCoeff(N,n_out,n_in); %Obtain boundary coefficients for SPN
        
        
        Cgrad=[(1+B(1))./(3*mua1(belem)), -D(1)./mua3(belem),-F(1)./mua5(belem),-H(1)./mua7(belem)
            -D(2)./mua1(belem), (1+B(2))./(7*mua3(belem)),-F(2)./mua5(belem),-H(2)./mua7(belem)
            -D(3)./mua1(belem),-F(3)./mua3(belem),(1+B(3))./(11*mua5(belem)),-H(3)./mua7(belem)
            -D(4)./mua1(belem),-F(4)./mua3(belem),-H(4)./mua5(belem),(1+B(4))./(15*mua7(belem))];
        
        C=[(0.5+A(1))*ones(1,1,nb), -(1/8+C(1))*ones(1,1,nb),-(-1/16+E(1))*ones(1,1,nb),-(5/128+G(1))*ones(1,1,nb)
            -(1/8+C(2))*ones(1,1,nb),(7/24+A(2))*ones(1,1,nb),-(41/384+E(2))*ones(1,1,nb),-(-1/16+G(2))*ones(1,1,nb)
            -(-1/16+C(3))*ones(1,1,nb),-(41/384+E(3))*ones(1,1,nb),(407/1920+A(3))*ones(1,1,nb),-(233/2560+G(3))*ones(1,1,nb)
            -(5/128+C(4))*ones(1,1,nb),-(-1/16+E(4))*ones(1,1,nb),-(233/2560+G(4))*ones(1,1,nb),(3023/17920+A(4))*ones(1,1,nb)];
        
        
end
end