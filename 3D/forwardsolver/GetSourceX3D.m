function S=GetSourceX3D(N,P,T,Src,B)
% Modified on 5th December 2020 to support extended sources
%Compute the source terms in the SP-3 formulation according to eqn (58b)
D = size(T,2);
Nodes = size(P,1);
NofSrcElems = size(Src.elem,1);
S=zeros(Nodes*(N+1)/2,NofSrcElems);


F1=@(y) 2*abs(cos(y)).*sin(y);
F2=@(x) (5*(abs(cos(x)).^3)-3*abs(cos(x))).*sin(x);
CS(1,1,1:NofSrcElems)=integral(F1,pi/2,pi)*Src.Strength*ones(1,1,NofSrcElems);  %check this integral. shouldnt it go from pi to pi/2 what difference does it make?
CS(1,2,1:NofSrcElems)=integral(F2,pi/2,pi)*Src.Strength*ones(1,1,NofSrcElems);

if N>3
    F3=@(x1) ((63/4)*abs(cos(x1)).^5 -(35/2)*abs(cos(x1)).^3+(15/4)*abs(cos(x1))).*sin(x1);
    CS(1,3,1:NofSrcElems)=integral(F3,pi/2,pi)*Src.Strength*ones(1,1,NofSrcElems);
    if N>5
        F4=@(x2) ((429/8)*abs(cos(x2)).^7 - (693/8)*abs(cos(x2)).^5 + (315/8)*abs(cos(x2)).^3-(35/8)*abs(cos(x2))).*sin(x2);
        CS(1,4,1:NofSrcElems)=integral(F4,pi/2,pi)*Src.Strength*ones(1,1,NofSrcElems);
    end
end

CS=repmat(CS,[(N+1)/2,1,1]);
Temp=sum(B.*CS,2);
Temp=repmat(Temp,[1,D,1]);
Nmat = ones(NofSrcElems,D);
Nmat = Nmat.*ismember(T(Src.elem,:),Src.node);

Ntemp(1,:,:)=(Nmat./repmat(sum(Nmat,2),1,D)).';
Ntemp=repmat(Ntemp,[(N+1)/2,1,1]);
Temp=Ntemp.*Temp;

for i=1:NofSrcElems
    S1=zeros(Nodes,1); 
    S1(T(Src.elem(i),:))=Temp(1,:,i);
    
    if N==1
        S = [S,S1]; 
    else
        S2=zeros(Nodes,1);
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
end
S = sum(S,2);


%     T=[DT(:,1),DT(:,2),DT(:,3),DT(:,4)];
% Nodes=size(P,1);
% %     Elem=size(T,1);
% %     [~,~,e_n,~]=GetEdgeLength(P,T);
% %     [~,e_n]=boundedges_element(P,T);
% S=[];
% Sx=Src.Loc(:,1); Sy=Src.Loc(:,2); Sz=Src.Loc(:,3);
% NofSources=size(Src.Loc,1);
% SInd=zeros(NofSources,1);Q=1;
% for i=1:NofSources
%     SInd(i)=find((edgeElem(:,1)==Src.elem(i)));  %Index
% end
% 
% x1=P(T(Src.elem,1),1); x2=P(T(Src.elem,2),1); x3=P(T(Src.elem,3),1); x4=P(T(Src.elem,4),1);
% y1=P(T(Src.elem,1),2); y2=P(T(Src.elem,2),2); y3=P(T(Src.elem,3),2); y4=P(T(Src.elem,4),2);
% z1=P(T(Src.elem,1),3); z2=P(T(Src.elem,2),3); z3=P(T(Src.elem,3),3); z4=P(T(Src.elem,4),3);
% 
% a1=x2.*(y3.*z4-y4.*z3) - x3.*(y2.*z4-y4.*z2) + x4.*(y2.*z3-y3.*z2);
% a2=-(x1.*(y3.*z4-y4.*z3) - x3.*(y1.*z4-y4.*z1) + x4.*(y1.*z3-y3.*z1));
% a3=x1.*(y2.*z4-y4.*z2) - x2.*(y1.*z4-y4.*z1) + x4.*(y1.*z2-y2.*z1);
% a4=-(x1.*(y2.*z3-y3.*z2) - x2.*(y1.*z3-y3.*z1) + x3.*(y1.*z2-y2.*z1));
% 
% b1=-((y3.*z4-y4.*z3) - (y2.*z4-y4.*z2) + (y2.*z3-y3.*z2));
% b2=(y3.*z4-y4.*z3) - (y1.*z4-y4.*z1) + (y1.*z3-y3.*z1);
% b3=-((y2.*z4-y4.*z2) - (y1.*z4-y4.*z1) + (y1.*z2-y2.*z1));
% b4=(y2.*z3-y3.*z2) - (y1.*z3-y3.*z1) + (y1.*z2-y2.*z1);
% 
% c1=(x3.*z4-x4.*z3) - (x2.*z4-x4.*z2) + (x2.*z3-x3.*z2);
% c2=-((x3.*z4-x4.*z3) - (x1.*z4-x4.*z1) + (x1.*z3-x3.*z1));
% c3=(x2.*z4-x4.*z2) - (x1.*z4-x4.*z1) + (x1.*z2-x2.*z1);
% c4=-((x2.*z3-x3.*z2) - (x1.*z3-x3.*z1) + (x1.*z2-x2.*z1));
% 
% d1=(y3.*x4-y4.*x3) - (y2.*x4-y4.*x2) + (y2.*x3-y3.*x2);
% d2=-((y3.*x4-y4.*x3) - (y1.*x4-y4.*x1) + (y1.*x3-y3.*x1));
% d3=(y2.*x4-y4.*x2) - (y1.*x4-y4.*x1) + (y1.*x2-y2.*x1);
% d4=-((y2.*x3-y3.*x2) - (y1.*x3-y3.*x1) + (y1.*x2-y2.*x1));
% V=zeros(NofSources,1);
% for s=1:NofSources
%     V(s)=(1/6)*det([1,1,1,1; x1(s),x2(s),x3(s),x4(s); y1(s),y2(s),y3(s),y4(s);z1(s),z2(s),z3(s),z4(s)]);
% end
% 
% N1=(a1+b1.*Sx+c1.*Sy+d1.*Sz)./(6*V);  N2=(a2+b2.*Sx+c2.*Sy+d2.*Sz)./(6*V);  N3=(a3+b3.*Sx+c3.*Sy+d3.*Sz)./(6*V); N4=(a4+b4.*Sx+c4.*Sy+d4.*Sz)./(6*V);
% 
% Mask=ones(NofSources,4);
% for i=1:NofSources
%     if edgeType(SInd(i))==1
%         Mask(i,4)=0;
%     elseif edgeType(SInd(i))==2
%         Mask(i,1)=0;
%     elseif edgeType(SInd(i))==3
%         Mask(i,3)=0;
%     else
%         Mask(i,2)=0;
%     end
% end
% 
% CS = zeros(1,(N+1)/2,NofSources);
% for s = 1:NofSources
% if Src(s).AOI==90
%     F1=@(y) 2*abs(cos(y)).*sin(y);
%     F2=@(x) (5*(abs(cos(x)).^3)-3*abs(cos(x))).*sin(x);
%     %     Q=zeros(Elem,1);
%     %     Q(src)=1; %Specifying source strength
%     
%     %     CS1=integral(F1,pi/2,pi)*Q;  %check this integral. shouldnt it go from pi to pi/2 what difference does it make?
%     %     CS2=integral(F2,pi/2,pi)*Q;
%     %
%     
%     CS(1,1,s)=2*pi*integral(F1,pi/2,pi)*Q;  %check this integral. shouldnt it go from pi to pi/2 what difference does it make?
%     CS(1,2,s)=2*pi*integral(F2,pi/2,pi)*Q;
%     
%     if N>3
%         F3=@(x1) ((63/4)*abs(cos(x1)).^5 -(35/2)*abs(cos(x1)).^3+(15/4)*abs(cos(x1))).*sin(x1);
%         CS(1,3,s)=2*pi*integral(F3,pi/2,pi)*Q;
%         if N>5
%             F4=@(x2) ((429/8)*abs(cos(x2)).^7 - (693/8)*abs(cos(x2)).^5 + (315/8)*abs(cos(x2)).^3-(35/8)*abs(cos(x2))).*sin(x2);
%             CS(1,4,s)=2*pi*integral(F4,pi/2,pi)*Q;
%         end
%     end
%     
% else
%     AOI=Src(s).AOI*pi/180;
% %     CS(1)=2*pi*(2*abs(cos(AOI)).*sin(AOI));
% % IMPT: we do not consider absolute value of the cosine here. This is
% % because we are considering the angle w.r.t to the outward normal. When
% % integrating over a set of angles however the abs() must be considered,
% % since then the integration limits are specified to allow only incoming
% % source directions. 
%       CS(1,1,s)=Src.Strength(s)*2*pi*(2*(cos(AOI)));
% %     CS(2)=2*pi*(5*(abs(cos(AOI)).^3)-3*abs(cos(AOI))).*sin(AOI);
%     CS(1,2,s)=Src.Strength(s)*2*pi*(5*((cos(AOI)).^3)-3*(cos(AOI)));
%     if N>3
%         CS(1,3,s)=Src.Strength(s)*2*pi*((63/4)*(cos(AOI)).^5 -(35/2)*(cos(AOI)).^3+(15/4)*(cos(AOI))).*sin(AOI);
%         if N>5
%             CS(1,4,s)=Src.Strength(s)*2*pi*((429/8)*(cos(AOI)).^7 - (693/8)*(cos(AOI)).^5 + (315/8)*(cos(AOI)).^3-(35/8)*(cos(AOI))).*sin(AOI);
%         end
%     end
% end
% 
% CS=repmat(CS,[(N+1)/2,1,1]);
% %     Temp=B(:,:,src)*[CS1(src);CS2(src)];  %Have changed this coz B is a (2x2xelement) matrix here
% Temp=sum(B.*CS,2);
% % Temp=CS';
% Temp=repmat(Temp,[1,4,1]);
% Ntemp=zeros(1,4,NofSources);
% Ntemp(1,:,:)=(Mask.*[N1,N2,N3,N4])';
% Ntemp=repmat(Ntemp,[(N+1)/2,1,1]);
% Temp=Ntemp.*Temp;
% for i=1:NofSources
%     S1=zeros(Nodes,1); S2=zeros(Nodes,1);
%     S1(T(Src.elem(i),:))=Temp(1,:,i);
%     S2(T(Src.elem(i),:))=Temp(2,:,i);
%     if N==3
%         S=[S,[S1;S2]];
%     else
%         S3=zeros(Nodes,1);
%         S3(T(Src.elem(i),:))=Temp(3,:,i);
%         if N==5
%             S=[S,[S1;S2;S3]];
%         else
%             S4=zeros(Nodes,1);
%             S4(T(Src.elem(i),:))=Temp(4,:,i);
%             S=[S,[S1;S2;S3;S4]];
%         end
%     end
% end
% 
% end