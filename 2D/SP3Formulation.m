function [Cgradphi,Cphi,Cgradbphi,Cbphi,CJ]=SP3Formulation(N,Nelem,mesh,mua,mus,bn,n0,ni,Src,Det)
%%
% function to compute the elemental coeffecient matrices for the SPN
% approximation. Presently only works for the SP3 approximation. Support
% can be extended to higher/ lower orders in future versions. (Not tested
% for n=1,5,7)

%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%
% N - order of the approximation
% Nelem - # of elements
% mesh - structure variable containing meshing information
% mua - elementwise total absorption coeffecient. for freq domain case, this includes
% the freq. dependednt term as well. 
% mus - elementwise scattering coeffecient
% g - elementwise anisotropy factor
% n0 - refractive index of the outside medium 
% ni - refractive index of the medium
% Src - structure containing the source information
% Det - structure containing the detector information

%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%
% Cgradphi,Cphi,Cgradbphi,Cbphi, CJ - see ref. 1 for definitions and
% details

%%%%%%%%% REFERENCES %%%%%%%%%%%
% 1. Fully non-linear Sp3 approximation based FOT - self
% 2. Klose 2005
%%%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
edgeElem = mesh.edgeElem;
BElem=size(edgeElem,1);
NofDet=numel(Det.elem);
NofSources = size(Src,2);
%Moments of absorption ref. Klose 2005

mua1=((mus*(1-bn(1)))+mua);
mua2=((mus*(1-bn(2)))+mua);
mua3=((mus*(1-bn(3)))+mua);
mua0=mua;
% Initialise matrices
Cgradbphi=zeros((N+1)/2,(N+1)/2,Nelem);
Cbphi=zeros((N+1)/2,(N+1)/2,Nelem);
Cj=zeros(1,(N+1)/2,Nelem);
Cjgrad=zeros(1,(N+1)/2,Nelem);
if N==1
    % This is DA and not SP1. Difference is primarily in boundary
    % condition. 
    [A,B]=CalcReflCoeff(N,n0,ni);
    Cgradphi=1./(3*mua1);
    Cphi = mua0;
    Cgradbphi(:,:,edgeElem(:,1))=(1+B(1))*Cgradphi(:,:,edgeElem(:,1));
    Cbphi(:,:,edgeElem(:,1)) = (0.5+A(1))*ones(1,1,BElem);
    Cj(:,:,Det.elem)=1;
    
else if N==3
        Cgradphi=[1./(3*mua1), zeros(1,1,Nelem); zeros(1,1,Nelem), 1./(7*mua3)];
        Cphi=[mua0, -(2/3)*mua0; -(2/3)*mua0, (4/9)*mua0+(5/9)*mua2];
        
        %Boundary matrices
        % [~,e_n]=boundedges_element(P,T);
        if length(n0)==1
            % if same refractive index on all external interfaces
        [A,B,C,D]=CalcReflCoeff(N,n0,ni); %Obtain boundary coefficients for SPN
        Cgradbphi(:,:,edgeElem(:,1))=[(1+B(1))./(3*mua1(edgeElem(:,1))), -D(1)./mua3(edgeElem(:,1))
            -D(2)./mua1(edgeElem(:,1)), (1+B(2))./(7*mua3(edgeElem(:,1)))];
        
        Cbphi(:,:,edgeElem(:,1))=[(0.5+A(1))*ones(1,1,BElem), -(1/8+C(1))*ones(1,1,BElem)
            -(1/8+C(2))*ones(1,1,BElem),(7/24+A(2))*ones(1,1,BElem)];
        else
            % if refractive index is different across different external interfaces
            [Cgradbphi,Cbphi]=DefineBoundMat(N,n0,ni,mua1,mua3,mesh);
        end
        

        % For boundary sources, fiber may have different refractive index, so
        % coefficients  may be different here
        for s=1:NofSources
            if (Src(s).ID~='I')
                NofSrc_elem = size(Src(s).elem,1);
                [A,B,C,D]=CalcReflCoeff(N,Src(s).nfibre,ni); %Obtain boundary coefficients for SPN
                Cgradbphi(:,:,Src(s).elem(:,1))=[(1+B(1))./(3*mua1(Src(s).elem(:,1))), -D(1)./mua3(Src(s).elem(:,1))
                    -D(2)./mua1(Src(s).elem(:,1)), (1+B(2))./(7*mua3(Src(s).elem(:,1)))];
                
                Cbphi(:,:,Src(s).elem(:,1))=[(0.5+A(1))*ones(1,1,NofSrc_elem), -(1/8+C(1))*ones(1,1,NofSrc_elem)
                    -(1/8+C(2))*ones(1,1,NofSrc_elem),(7/24+A(2))*ones(1,1,NofSrc_elem)];
            end
        end
        %Matrices for partial current evaluation
        %Note: these are of size 1x2 since the partial current is expressed as a linear
        %combination of the composite moments
        
        if Det.AOD==90
            J=CalcDetCoeff(N,Det.nfibre,ni);
            Cj(:,:,Det.elem)=[((1/4)+J(1))*ones(1,1,NofDet), (-(2/3)*((1/4)+J(1))+(1/3)*((5/16)+J(3)))*ones(1,1,NofDet)];
            Cjgrad(:,:,Det.elem)=[(0.5+J(2))./(3*mua1(Det.elem)), J(4)./(7*mua3(Det.elem))];
        else
            J=DirectedDetection(N,Det.nfibre,ni,Det.AOD);
            Cj(:,:,Det.elem)=[(J(1))*ones(1,1,NofDet), J(3)*ones(1,1,NofDet)];
            Cjgrad(:,:,Det.elem)=[J(2)./(2*mua1(Det.elem)), J(4)./(2*mua3(Det.elem))];
        end
        
    else if N==5
%             mua4=((mus.*(1-g.^4))+mua);
%             mua5=((mus.*(1-g.^5))+mua);
%             Cgradphi=[1./(3*mua1), zeros(1,1,Nelem), zeros(1,1,Nelem)
%                 zeros(1,1,Nelem), 1./(7*mua3),zeros(1,1,Nelem)
%                 zeros(1,1,Nelem), zeros(1,1,Nelem),1./(11*mua5)];
%             
%             Cphi=[  mua0, -(2/3)*mua0, (8/15)*mua0
%                 -(2/3)*mua0, (4/9)*mua0+(5/9)*mua2, -((16/45)*mua0 + (4/9)*mua2)
%                 (8/15)*mua0,-((16/45)*mua0 + (4/9)*mua2), ((64/225)*mua0 + (16/45)*mua2 + (9/25)*mua4)];
%             
%             %Boundary Matrices
%             [A,B,C,D,E,F]=CalcReflCoeff(N,n0,ni); %Obtain boundary coefficients for SPN
%             
%             
%             Cgradbphi(:,:,edgeElem(:,1))=[(1+B(1))./(3*mua1(edgeElem(:,1))), -D(1)./mua3(edgeElem(:,1)),-F(1)./mua5(edgeElem(:,1))
%                 -D(2)./mua1(edgeElem(:,1)), (1+B(2))./(7*mua3(edgeElem(:,1))),-F(2)./mua5(edgeElem(:,1))
%                 -D(3)./mua1(edgeElem(:,1)),-F(3)./mua3(edgeElem(:,1)),(1+B(3))./(11*mua5(edgeElem(:,1)))];
%             
%             Cbphi(:,:,edgeElem(:,1))=[(0.5+A(1))*ones(1,1,BElem), -(1/8+C(1))*ones(1,1,BElem),-(-1/16+E(1))*ones(1,1,BElem)
%                 -(1/8+C(2))*ones(1,1,BElem),(7/24+A(2))*ones(1,1,BElem),-(41/384+E(2))*ones(1,1,BElem)
%                 -(-1/16+C(3))*ones(1,1,BElem),-(41/384+E(3))*ones(1,1,BElem),(407/1920+A(3))*ones(1,1,BElem)];
%             
%             % Value at Source elem location
%             if (Src.ID~='I')
%             [A,B,C,D,E,F]=CalcReflCoeff(N,Src.nfibre,ni); %Obtain boundary coefficients for SPN
%             
%             
%             Cgradbphi(:,:,Src.elem(:,1))=[(1+B(1))./(3*mua1(Src.elem(:,1))), -D(1)./mua3(Src.elem(:,1)),-F(1)./mua5(Src.elem(:,1))
%                 -D(2)./mua1(Src.elem(:,1)), (1+B(2))./(7*mua3(Src.elem(:,1))),-F(2)./mua5(Src.elem(:,1))
%                 -D(3)./mua1(Src.elem(:,1)),-F(3)./mua3(Src.elem(:,1)),(1+B(3))./(11*mua5(Src.elem(:,1)))];
%             
%             Cbphi(:,:,Src.elem(:,1))=[(0.5+A(1))*ones(1,1,NofSources), -(1/8+C(1))*ones(1,1,NofSources),-(-1/16+E(1))*ones(1,1,NofSources)
%                 -(1/8+C(2))*ones(1,1,NofSources),(7/24+A(2))*ones(1,1,NofSources),-(41/384+E(2))*ones(1,1,NofSources)
%                 -(-1/16+C(3))*ones(1,1,NofSources),-(41/384+E(3))*ones(1,1,NofSources),(407/1920+A(3))*ones(1,1,NofSources)];
%             end
%             
%             %Partial Current Matrices
%             if Det.AOD==90
%                 J=CalcDetCoeff(N,Det.nfibre,ni);
%                 Cj(:,:,Det.elem)=[((1/4)+J(1))*ones(1,1,NofDet), (-(2/3)*((1/4)+J(1))+(1/3)*((5/16)+J(3)))*ones(1,1,NofDet),((8/15)*((1/4)+J(1))-(4/15)*((5/16)+J(3))+(1/5)*(-(3/32)+J(5)))*ones(1,1,NofDet)];
%                 Cjgrad(:,:,Det.elem)=[(0.5+J(2))./(3*mua1(Det.elem)), J(4)./(7*mua3(Det.elem)), J(6)./(11*mua5(Det.elem))];
%             else
%                 J=DirectedDetection(N,Det.nfibre,ni,Det.AOD);
%                 Cj(:,:,Det.elem)=[(J(1))*ones(1,1,NofDet), J(3)*ones(1,1,NofDet),J(5)*ones(1,1,NofDet)];
%                 Cjgrad(:,:,Det.elem)=[J(2)./(2*mua1(Det.elem)), J(4)./(2*mua3(Det.elem)),J(6)./(2*mua5(Det.elem))];
%             end
        else
            %(N==7)
%             mua4=((mus.*(1-g.^4))+mua);
%             mua5=((mus.*(1-g.^5))+mua);
%             mua6=((mus.*(1-g.^6))+mua);
%             mua7=((mus.*(1-g.^7))+mua);
%             Cgradphi=[1./(3*mua1), zeros(1,1,Nelem), zeros(1,1,Nelem), zeros(1,1,Nelem)
%                 zeros(1,1,Nelem), 1./(7*mua3),zeros(1,1,Nelem),zeros(1,1,Nelem)
%                 zeros(1,1,Nelem), zeros(1,1,Nelem),1./(11*mua5),zeros(1,1,Nelem)
%                 zeros(1,1,Nelem),zeros(1,1,Nelem),zeros(1,1,Nelem),1./(15*mua7)];
%             Cphi=[  mua0, -(2/3)*mua0, (8/15)*mua0, -(16/35)*mua0
%                 -(2/3)*mua0, (4/9)*mua0+(5/9)*mua2, -((16/45)*mua0 + (4/9)*mua2), ((32/105)*mua0 + (8/21)*mua2)
%                 (8/15)*mua0,-((16/45)*mua0 + (4/9)*mua2), ((64/225)*mua0 + (16/45)*mua2 + (9/25)*mua4), -((128/525)*mua0 + (32/105)*mua2 + (54/175)*mua4)
%                 -(16/35)*mua0,((32/105)*mua0 + (8/21)*mua2),-((128/525)*mua0 + (32/105)*mua2 + (54/175)*mua4), ((256/1225)*mua0+(64/245)*mua2 + (324/1225)*mua4 + (13/49)*mua6)];
%             
%             %Boundary matrices
%             [A,B,C,D,E,F,G,H]=CalcReflCoeff(N,n0,ni); %Obtain boundary coefficients for SPN
%             
%             
%             Cgradbphi(:,:,edgeElem(:,1))=[(1+B(1))./(3*mua1(edgeElem(:,1))), -D(1)./mua3(edgeElem(:,1)),-F(1)./mua5(edgeElem(:,1)),-H(1)./mua7(edgeElem(:,1))
%                 -D(2)./mua1(edgeElem(:,1)), (1+B(2))./(7*mua3(edgeElem(:,1))),-F(2)./mua5(edgeElem(:,1)),-H(2)./mua7(edgeElem(:,1))
%                 -D(3)./mua1(edgeElem(:,1)),-F(3)./mua3(edgeElem(:,1)),(1+B(3))./(11*mua5(edgeElem(:,1))),-H(3)./mua7(edgeElem(:,1))
%                 -D(4)./mua1(edgeElem(:,1)),-F(4)./mua3(edgeElem(:,1)),-H(4)./mua5(edgeElem(:,1)),(1+B(4))./(15*mua7(edgeElem(:,1)))];
%             
%             Cbphi(:,:,edgeElem(:,1))=[(0.5+A(1))*ones(1,1,BElem), -(1/8+C(1))*ones(1,1,BElem),-(-1/16+E(1))*ones(1,1,BElem),-(5/128+G(1))*ones(1,1,BElem)
%                 -(1/8+C(2))*ones(1,1,BElem),(7/24+A(2))*ones(1,1,BElem),-(41/384+E(2))*ones(1,1,BElem),-(-1/16+G(2))*ones(1,1,BElem)
%                 -(-1/16+C(3))*ones(1,1,BElem),-(41/384+E(3))*ones(1,1,BElem),(407/1920+A(3))*ones(1,1,BElem),-(233/2560+G(3))*ones(1,1,BElem)
%                 -(5/128+C(4))*ones(1,1,BElem),-(-1/16+E(4))*ones(1,1,BElem),-(233/2560+G(4))*ones(1,1,BElem),(3023/17920+A(4))*ones(1,1,BElem)];
%             
%             % At source location
%              [A,B,C,D,E,F,G,H]=CalcReflCoeff(N,Src.nfibre,ni); %Obtain boundary coefficients for SPN
%             
%             
%             Cgradbphi(:,:,Src.elem(:,1))=[(1+B(1))./(3*mua1(Src.elem(:,1))), -D(1)./mua3(Src.elem(:,1)),-F(1)./mua5(Src.elem(:,1)),-H(1)./mua7(Src.elem(:,1))
%                 -D(2)./mua1(Src.elem(:,1)), (1+B(2))./(7*mua3(Src.elem(:,1))),-F(2)./mua5(Src.elem(:,1)),-H(2)./mua7(Src.elem(:,1))
%                 -D(3)./mua1(Src.elem(:,1)),-F(3)./mua3(Src.elem(:,1)),(1+B(3))./(11*mua5(Src.elem(:,1))),-H(3)./mua7(Src.elem(:,1))
%                 -D(4)./mua1(Src.elem(:,1)),-F(4)./mua3(Src.elem(:,1)),-H(4)./mua5(Src.elem(:,1)),(1+B(4))./(15*mua7(Src.elem(:,1)))];
%             
%             Cbphi(:,:,Src.elem(:,1))=[(0.5+A(1))*ones(1,1,NofSources), -(1/8+C(1))*ones(1,1,NofSources),-(-1/16+E(1))*ones(1,1,NofSources),-(5/128+G(1))*ones(1,1,NofSources)
%                 -(1/8+C(2))*ones(1,1,NofSources),(7/24+A(2))*ones(1,1,NofSources),-(41/384+E(2))*ones(1,1,NofSources),-(-1/16+G(2))*ones(1,1,NofSources)
%                 -(-1/16+C(3))*ones(1,1,NofSources),-(41/384+E(3))*ones(1,1,NofSources),(407/1920+A(3))*ones(1,1,NofSources),-(233/2560+G(3))*ones(1,1,NofSources)
%                 -(5/128+C(4))*ones(1,1,NofSources),-(-1/16+E(4))*ones(1,1,NofSources),-(233/2560+G(4))*ones(1,1,NofSources),(3023/17920+A(4))*ones(1,1,NofSources)];
%             
%             
%             
%             %Partial Current matrices
%             if Det.AOD==90
%                 J=CalcDetCoeff(N,Det.nfibre,ni);
%                 Cj(:,:,Det.elem)=[((1/4)+J(1))*ones(1,1,NofDet), (-(2/3)*((1/4)+J(1))+(1/3)*((5/16)+J(3)))*ones(1,1,NofDet),((8/15)*((1/4)+J(1))-(4/15)*((5/16)+J(3))+(1/5)*(-(3/32)+J(5)))*ones(1,1,NofDet),(-(16/35)*((1/4)+J(1))+(8/35)*((5/16)+J(3))-(6/35)*(-(3/32)+J(5))+(1/7)*((13/256+J(7))))*ones(1,1,NofDet)];
%                 Cjgrad(:,:,Det.elem)=[(0.5+J(2))./(3*mua1(Det.elem)), J(4)./(7*mua3(Det.elem)),J(6)./(11*mua5(Det.elem)),J(8)/(15*mua7(Det.elem))];
%             else
%                 J=DirectedDetection(N,Det.nfibre,ni,Det.AOD);
%                 Cj(:,:,Det.elem)=[(J(1))*ones(1,1,NofDet), J(3)*ones(1,1,NofDet),J(5)*ones(1,1,NofDet),J(7)*ones(1,1,NofDet)];
%                 Cjgrad(:,:,Det.elem)=[J(2)./(2*mua1(Det.elem)), J(4)./(2*mua3(Det.elem)),J(6)./(2*mua5(Det.elem)),J(8)./(2*mua7(Det.elem))];
%             end
%             
        end
        
    end
end

CJ=zeros(1,(N+1)/2,NofDet);
for i=1:NofDet
    Temp=Cgradbphi(:,:,Det.elem(i))\Cbphi(:,:,Det.elem(i));
    CJ(:,:,i)=Cj(:,:,Det.elem(i))+Cjgrad(:,:,Det.elem(i))*Temp;
end

end
