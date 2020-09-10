function S=GetSourceX2D(N,P,T,edgeElem,edgeType,Src,B)
% Created on 9th April 2019
Nodes=size(P,1);
S=[];
Sx=Src.Loc(:,1); Sy=Src.Loc(:,2);
NofSources=size(Src.Loc,1);
SInd=zeros(NofSources,1);
for i=1:NofSources
    SInd(i)=find((edgeElem(:,1)==Src.elem(i)));  %Index
end

x1=P(T(Src.elem,1),1); x2=P(T(Src.elem,2),1); x3=P(T(Src.elem,3),1);
y1=P(T(Src.elem,1),2); y2=P(T(Src.elem,2),2); y3=P(T(Src.elem,3),2);

a1=x2.*y3-y2.*x3;   a2=x3.*y1-y3.*x1;   a3=x1.*y2-y1.*x2;
b1=y2-y3;           b2=y3-y1;           b3=y1-y2;
c1=x3-x2;           c2=x1-x3;           c3=x2-x1;
A=0.5*(b1.*c2-b2.*c1);

N1=(a1+b1.*Sx+c1.*Sy)./(2*A);  N2=(a2+b2.*Sx+c2.*Sy)./(2*A);  N3=(a3+b3.*Sx+c3.*Sy)./(2*A);

Mask=ones(NofSources,3);
for i=1:NofSources
    if edgeType(SInd(i))==1
        Mask(i,3)=0;
%                Mask(i,1:2)=sqrt((x1(i)-x2(i)).^2 + (y1(i)-y2(i)).^2);
    else if edgeType(SInd(i))==2
            Mask(i,1)=0;
%                        Mask(i,2:3)=sqrt((x3(i)-x2(i)).^2 + (y3(i)-y2(i)).^2);
        else
            Mask(i,2)=0;
%                        Mask(i,[1,3])=sqrt((x1(i)-x3(i)).^2 + (y1(i)-y3(i)).^2);
        end
    end
end
CS = zeros(1,(N+1)/2,NofSources);
for s = 1:NofSources
    if Src.AOI==90 %Isotropic source.
        
        F1=@(y) 2*abs(cos(y)).*sin(y);
        F2=@(x) (5*(abs(cos(x)).^3)-3*abs(cos(x))).*sin(x);
        CS(1,1,s)=integral(F1,pi/2,pi)*Src.Strength(s);  %check this integral. shouldnt it go from pi to pi/2 what difference does it make?
        CS(1,2,s)=integral(F2,pi/2,pi)*Src.Strength(s);
        
        if N>3
            F3=@(x1) ((63/4)*abs(cos(x1)).^5 -(35/2)*abs(cos(x1)).^3+(15/4)*abs(cos(x1))).*sin(x1);
            CS(1,3,s)=integral(F3,pi/2,pi)*Src.Strength(s);
            if N>5
                F4=@(x2) ((429/8)*abs(cos(x2)).^7 - (693/8)*abs(cos(x2)).^5 + (315/8)*abs(cos(x2)).^3-(35/8)*abs(cos(x2))).*sin(x2);
                CS(1,4,s)=integral(F4,pi/2,pi)*Src.Strength(s);
            end
        end
        
    else
        
        AOI=Src.AOI*pi/180;
%         if Src.NA ==0
            CS(1,1,s)=Src.Strength(s)*2*abs(cos(AOI)); % Check definition of delta function as in Ishimaru. Ch 7 Pg 160
            CS(1,2,s)=Src.Strength(s)*(5*(abs(cos(AOI)).^3)-3*abs(cos(AOI)));
            if N>3
                CS(1,3,s)=Src.Strength(s)*((63/4)*abs(cos(AOI)).^5 -(35/2)*abs(cos(AOI)).^3+(15/4)*abs(cos(AOI)));
                if N>5
                    CS(1,4,s)=Src.Strength(s)*((429/8)*abs(cos(AOI)).^7 - (693/8)*abs(cos(AOI)).^5 + (315/8)*abs(cos(AOI)).^3-(35/8)*abs(cos(AOI)));
                end
            end
%         else
%             F1=@(y) 2*abs(cos(y)).*sin(y);
%         F2=@(x) (5*(abs(cos(x)).^3)-3*abs(cos(x))).*sin(x);
%         CS(1,1,s)=integral(F1,AOI-Src.NA/2,AOI+Src.NA/2)*Src.Strength(s);  %check this integral. shouldnt it go from pi to pi/2 what difference does it make?
%         CS(1,2,s)=integral(F2,AOI-Src.NA/2,AOI+Src.NA/2)*Src.Strength(s);
%         
%         if N>3
%             F3=@(x1) ((63/4)*abs(cos(x1)).^5 -(35/2)*abs(cos(x1)).^3+(15/4)*abs(cos(x1))).*sin(x1);
%             CS(1,3,s)=integral(F3,pi/2,pi)*Src.Strength(s);
%             if N>5
%                 F4=@(x2) ((429/8)*abs(cos(x2)).^7 - (693/8)*abs(cos(x2)).^5 + (315/8)*abs(cos(x2)).^3-(35/8)*abs(cos(x2))).*sin(x2);
%                 CS(1,4,s)=integral(F4,pi/2,pi)*Src.Strength(s);
%             end
%         end
%         end
        
    end
end

CS=repmat(CS,[(N+1)/2,1,1]);
Temp=sum(B.*CS,2);
Temp=repmat(Temp,[1,3,1]);
Ntemp=zeros(1,3,NofSources);
Ntemp(1,:,:)=(Mask.*[N1,N2,N3])';
Ntemp=repmat(Ntemp,[(N+1)/2,1,1]);
Temp=Ntemp.*Temp;
for i=1:NofSources
    S1=zeros(Nodes,1); S2=zeros(Nodes,1);
    S1(T(Src.elem(i),:))=Temp(1,:,i);
    S2(T(Src.elem(i),:))=Temp(2,:,i);
    if N==3
        S=[S,[S1;S2]];
    else
        S3=zeros(Nodes,1);
        S3(T(Src.elem(i),:))=Temp(3,:,i);
        if N==5
            S=[S,[S1;S2;S3]];
        else
            S4=zeros(Nodes,1);
            S4(T(Src.elem(i),:))=Temp(4,:,i);
            S=[S,[S1;S2;S3;S4]];
        end
    end
end


