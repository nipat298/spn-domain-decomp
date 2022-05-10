function [dCgrad,dC,dCgradB,dCjgrad]=dSP3Formulation(N,mesh,mua,mus,bn,n0,ni,Det,varargin)
%the derivative matrices for the SP-n matrices
% Though we compute dcjgrad, it is not used in the codes. If needed the
% section can be commented out to avoid confusion.

% NofElem = size(mesh.tri,1);
NofElem = length(mua);
edgeElem = mesh.edgeElem;
dmu=ones(1,1,NofElem);
% dmu(1,1,ENEL)=100;
if nargin>8
    dmu=varargin{1}*dmu;
end

dCgradB=zeros((N+1)/2,(N+1)/2,NofElem);
dCjgrad=zeros(1,(N+1)/2,NofElem);

mua1=((mus*(1-bn(2)))+mua);
mua3=((mus*(1-bn(4)))+mua);
if N==1
    
    dCgrad=-(1./(3*(mua1).^2)).*dmu;
    dC=ones(1,1,NofElem).*dmu;
    if length(n0)==1
        [~,B]=CalcReflCoeff(N,n0,ni);
        dCgradB(:,:,edgeElem)=-((1+B(1))./(3*(mua1(:,:,edgeElem)).^2)).*dmu(:,:,edgeElem);
    else
        dCgradB=Define_dBoundMat(N,n0,ni,mesh,dmu,mua1);
    end
    
elseif N==3
    
    dCgrad=[-(1./(3*(mua1.^2))).*dmu, zeros(1,1,NofElem); zeros(1,1,NofElem), -(1./(7*(mua3.^2))).*dmu];
    dC=[ones(1,1,NofElem).*dmu, -(2/3)*ones(1,1,NofElem).*dmu; -(2/3)*ones(1,1,NofElem).*dmu, ones(1,1,NofElem).*dmu];
    
    if length(n0)==1
        %Boundary matrix
        [~,B,~,D]=CalcReflCoeff(N,n0,ni);
        dCgradB(:,:,edgeElem)=[-((1+B(1))./(3*(mua1(edgeElem).^2))).*dmu(edgeElem), (D(1)./(mua3(edgeElem).^2)).*dmu(edgeElem); (D(2)./(mua1(edgeElem).^2)).*dmu(edgeElem), -((1+B(2))./(7*(mua3(edgeElem).^2))).*dmu(edgeElem)];
    else
        dCgradB=Define_dBoundMat(N,n0,ni,mesh,dmu,mua1,mua3);
    end
    %Partial Current matrix
    %         if Det.AOD==90
    %             J=CalcDetCoeff(N,Det.nfibre,ni);
    %             dCjgrad(:,:,Det.elem)=[-(0.5+J(2))./(3*((mua1(Det.elem).^2))), -J(4)./(7*((mua3(Det.elem).^2)))];
    %         else
    %             J=Directed.Detection(N,Det.nfibre,ni,Det.AOD);
    %             dCjgrad(:,:,Det.elem)=[-J(2)./(2*(mua1(Det.elem)).^2), -J(4)./(2*(mua3(Det.elem)).^2)];
    %         end
elseif N==5
    mua5=((mus.*(1-bn(6)))+mua);
    dCgrad=[-(1./(3*(mua1.^2))).*dmu, zeros(1,1,NofElem),zeros(1,1,NofElem)
        zeros(1,1,NofElem), -(1./(7*(mua3.^2))).*dmu,zeros(1,1,NofElem)
        zeros(1,1,NofElem),zeros(1,1,NofElem),-(1./(11*(mua5.^2))).*dmu];
    dC=[ones(1,1,NofElem).*dmu, -(2/3)*ones(1,1,NofElem).*dmu, (8/15)*ones(1,1,NofElem).*dmu
        -(2/3)*ones(1,1,NofElem).*dmu, ones(1,1,NofElem).*dmu, -(4/5)*ones(1,1,NofElem).*dmu
        (8/15)*ones(1,1,NofElem).*dmu,-(4/5)*ones(1,1,NofElem).*dmu, ones(1,1,NofElem).*dmu];
    %             %Boundary matrix
    if (length(n0)==1)
        [~,B,~,D,~,F]=CalcReflCoeff(N,n0,ni);
        dCgradB(:,:,edgeElem)=[-((1+B(1))./(3*(mua1(edgeElem).^2))).*dmu(edgeElem), (D(1)./(mua3(edgeElem).^2)).*dmu(edgeElem),(F(1)./mua5(edgeElem).^2).*dmu(edgeElem)
            (D(2)./(mua1(edgeElem).^2)).*dmu(edgeElem), -((1+B(2))./(7*(mua3(edgeElem).^2))).*dmu(edgeElem),(F(2)./mua5(edgeElem).^2).*dmu(edgeElem)
            (D(3)./(mua1(edgeElem).^2)).*dmu(edgeElem), (F(3)./(mua3(edgeElem).^2)).*dmu(edgeElem),-((1+B(3))./(11*(mua5(edgeElem).^2))).*dmu(edgeElem)];
    else
        dCgradB= Define_dBoundMat(N,n0,ni,mua1,mua3,mua5,mesh,dmu);
    end
    %             %Partial Current matrix
    %             if Det.AOD==90
    %                 J=CalcDetCoeff(N,Det.nfibre,ni);
    %                 dCjgrad(:,:,Det.elem)=[-(0.5+J(2))./(3*((mua1(Det.elem).^2))), -J(4)./(7*((mua3(Det.elem).^2))), -J(6)./(11.*(mua5(Det.elem).^2))];
    %             else
    %                 J=Directed.Detection(N,Det.nfibre,ni,Det.AOD);
    %                 dCjgrad(:,:,Det.elem)=[-J(2)./(2*(mua1(Det.elem)).^2), -J(4)./(2*(mua3(Det.elem)).^2), -J(6)./(2*(mua5(Det.elem)).^2)];
    %             end
else
    mua5=((mus.*(1-bn(6)))+mua);
    mua7=((mus.*(1-bn(8)))+mua);
    dCgrad=[-(1./(3*(mua1.^2))).*dmu, zeros(1,1,NofElem),zeros(1,1,NofElem),zeros(1,1,NofElem)
        zeros(1,1,NofElem), -(1./(7*(mua3.^2))).*dmu,zeros(1,1,NofElem),zeros(1,1,NofElem)
        zeros(1,1,NofElem),zeros(1,1,NofElem),-(1./(11*(mua5.^2))).*dmu,zeros(1,1,NofElem)
        zeros(1,1,NofElem),zeros(1,1,NofElem),zeros(1,1,NofElem),-(1./(15*(mua7.^2))).*dmu];
    dC=[ones(1,1,NofElem).*dmu, -(2/3)*ones(1,1,NofElem).*dmu, (8/15)*ones(1,1,NofElem).*dmu, -(16/35)*ones(1,1,NofElem).*dmu
        -(2/3)*ones(1,1,NofElem).*dmu, ones(1,1,NofElem).*dmu, -(4/5)*ones(1,1,NofElem).*dmu, (24/35)*ones(1,1,NofElem).*dmu
        (8/15)*ones(1,1,NofElem).*dmu,-(4/5)*ones(1,1,NofElem).*dmu, ones(1,1,NofElem).*dmu,-(30/29)*ones(1,1,NofElem).*dmu
        -(16/35)*ones(1,1,NofElem).*dmu,(24/35)*ones(1,1,NofElem).*dmu,-(30/29)*ones(1,1,NofElem).*dmu,ones(1,1,NofElem).*dmu];
    %             %Boundary matrix
    if (length(n0)==1)
        [~,B,~,D,~,F,~,H]=CalcReflCoeff(N,n0,ni);
        dCgradB(:,:,edgeElem)=[-((1+B(1))./(3*(mua1(edgeElem).^2))).*dmu(edgeElem), (D(1)./(mua3(edgeElem).^2)).*dmu(edgeElem),(F(1)./mua5(edgeElem).^2).*dmu(edgeElem),(H(1)./mua7(edgeElem).^2).*dmu(edgeElem)
            (D(2)./(mua1(edgeElem).^2)).*dmu(edgeElem), -((1+B(2))./(7*(mua3(edgeElem).^2))).*dmu(edgeElem),(F(2)./mua5(edgeElem).^2).*dmu(edgeElem),(H(2)./mua7(edgeElem).^2).*dmu(edgeElem)
            (D(3)./(mua1(edgeElem).^2)).*dmu(edgeElem), (F(3)./(mua3(edgeElem).^2)).*dmu(edgeElem),-((1+B(3))./(11*(mua5(edgeElem).^2))).*dmu(edgeElem),(H(3)./mua7(edgeElem).^2).*dmu(edgeElem)
            (D(4)./(mua1(edgeElem).^2)).*dmu(edgeElem), (F(4)./(mua3(edgeElem).^2)).*dmu(edgeElem),(H(4)./mua5(edgeElem).^2).*dmu(edgeElem),-((1+B(4))./(11*(mua7(edgeElem).^2))).*dmu(edgeElem)];
    else
        dCgradB=Define_dBoundMat(N,n0,ni,mua1,mua3,mua5,mua7,mesh,dmu);
    end
    %             %Partial Current matrix
    %             if Det.AOD==90
    %                 J=CalcDetCoeff(N,Det.nfibre,ni);
    %                 dCjgrad(:,:,Det.elem)=[-(0.5+J(2))./(3*((mua1(Det.elem).^2))), -J(4)./(7*((mua3(Det.elem).^2))), -J(6)./(11.*(mua5(Det.elem).^2)),-J(8)/(15*mua7(Det.elem).^2)];
    %             else
    %                 J=Directed.Detection(N,Det.nfibre,ni,Det.AOD);
    %                 dCjgrad(:,:,Det.elem)=[-J(2)./(2*(mua1(Det.elem)).^2), -J(4)./(2*(mua3(Det.elem)).^2), -J(6)./(2*(mua5(Det.elem)).^2),-J(8)./(2*(mua7(Det.elem)).^2)];
    %             end
end


end
