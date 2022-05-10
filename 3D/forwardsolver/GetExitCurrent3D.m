function RxFn=GetExitCurrent3D(N,CJ,P,T,Det)
% Modified on 5 Decemeber 2020 - supports extended detectors
%Compute the CJ term in equation (2)

NofDet=numel(Det.elem);
D = size(T,2);
%     T=[DT(:,1),DT(:,2),DT(:,3)];
Nodes=size(P,1);
%     Elem=size(T,1);

%     [~,~,e_n,~]=GetEdgeLength(P,T);
RxFn = zeros(NofDet,Nodes*(N+1)/2);
for i=1:NofDet
    if iscell(Det.elem)
        detelems = Det.elem{i};
        detnode = Det.node{i};
    else
        detelems = Det.elem(i);
        detnode = Det.node(i);
    end
    nelem = length(detelems);
    Nmat = ones(nelem,D);
    Nmat = Nmat.*ismember(T(detelems,:),detnode);
    Ntemp(1,:,:)=(Nmat./repmat(sum(Nmat,2),1,D)).';
    Ntemp=repmat(Ntemp,[(N+1)/2,1,1]);
    
    
    Temp=CJ(i,:,1:nelem);
    Temp = repmat(permute(Temp,[2,1,3]),1,D,1).*Ntemp;
    D1=zeros(Nodes,nelem);
    D2=D1; D3 = D1; D4 = D1;
    for j = 1:nelem
        switch N
            case 1
                D1(T(detelems(j),:),j)=Temp(1,:,j);
            case 3
                D1(T(detelems(j),:),j)=Temp(1,:,j);
                D2(T(detelems(j),:),j)=Temp(2,:,j);
            case 5
                D1(T(detelems(j),:),j)=Temp(1,:,j);
                D2(T(detelems(j),:),j)=Temp(2,:,j);
                D3(T(detelems(j),:),j)=Temp(3,:,j);
            case 7
                D1(T(detelems(j),:),j)=Temp(1,:,j);
                D2(T(detelems(j),:),j)=Temp(2,:,j);
                D3(T(detelems(j),:),j)=Temp(3,:,j);
                D4(T(detelems(j),:),j)=Temp(4,:,j);
        end
    end
    switch N
        
        case 1
            RxFn(i,:) =sum(D1,2).';
        case 3
            RxFn(i,:) =[sum(D1,2).',sum(D2,2).'];
        case 5
            RxFn(i,:) =[sum(D1,2).',sum(D2,2).',sum(D3,2).'];
        case 7
            RxFn(i,:) =[sum(D1,2).',sum(D2,2).',sum(D3,2).',sum(D4,2).'];
    end
    clear Ntemp Temp
end

end
% Dx=Det.Loc(:,1); Dy=Det.Loc(:,2);Dz=Det.Loc(:,3);
% %     T=[DT(:,1),DT(:,2),DT(:,3),DT(:,4)];
% Nodes=size(P,1);
% %     Elem=size(T,1);
% RxFn=zeros(d,((N+1)/2)*Nodes);
% %     [~,e_n]=boundedges_element(P,T);
% for i=1:d
%     DInd=edgeElem(:,1)==Det.elem(i);  %Index
%
%     x1=P(T(Det.elem(i),1),1); x2=P(T(Det.elem(i),2),1); x3=P(T(Det.elem(i),3),1); x4=P(T(Det.elem(i),4),1);
%     y1=P(T(Det.elem(i),1),2); y2=P(T(Det.elem(i),2),2); y3=P(T(Det.elem(i),3),2); y4=P(T(Det.elem(i),4),2);
%     z1=P(T(Det.elem(i),1),3); z2=P(T(Det.elem(i),2),3); z3=P(T(Det.elem(i),3),3); z4=P(T(Det.elem(i),4),3);
%
%     a1=x2.*(y3.*z4-y4.*z3) - x3.*(y2.*z4-y4.*z2) + x4.*(y2.*z3-y3.*z2);
%     a2=-(x1.*(y3.*z4-y4.*z3) - x3.*(y1.*z4-y4.*z1) + x4.*(y1.*z3-y3.*z1));
%     a3=x1.*(y2.*z4-y4.*z2) - x2.*(y1.*z4-y4.*z1) + x4.*(y1.*z2-y2.*z1);
%     a4=-(x1.*(y2.*z3-y3.*z2) - x2.*(y1.*z3-y3.*z1) + x3.*(y1.*z2-y2.*z1));
%
%     b1=-((y3.*z4-y4.*z3) - (y2.*z4-y4.*z2) + (y2.*z3-y3.*z2));
%     b2=(y3.*z4-y4.*z3) - (y1.*z4-y4.*z1) + (y1.*z3-y3.*z1);
%     b3=-((y2.*z4-y4.*z2) - (y1.*z4-y4.*z1) + (y1.*z2-y2.*z1));
%     b4=(y2.*z3-y3.*z2) - (y1.*z3-y3.*z1) + (y1.*z2-y2.*z1);
%
%     c1=(x3.*z4-x4.*z3) - (x2.*z4-x4.*z2) + (x2.*z3-x3.*z2);
%     c2=-((x3.*z4-x4.*z3) - (x1.*z4-x4.*z1) + (x1.*z3-x3.*z1));
%     c3=(x2.*z4-x4.*z2) - (x1.*z4-x4.*z1) + (x1.*z2-x2.*z1);
%     c4=-((x2.*z3-x3.*z2) - (x1.*z3-x3.*z1) + (x1.*z2-x2.*z1));
%
%     d1=(y3.*x4-y4.*x3) - (y2.*x4-y4.*x2) + (y2.*x3-y3.*x2);
%     d2=-((y3.*x4-y4.*x3) - (y1.*x4-y4.*x1) + (y1.*x3-y3.*x1));
%     d3=(y2.*x4-y4.*x2) - (y1.*x4-y4.*x1) + (y1.*x2-y2.*x1);
%     d4=-((y2.*x3-y3.*x2) - (y1.*x3-y3.*x1) + (y1.*x2-y2.*x1));
%
%     V=(1/6)*det([1,1,1,1; x1,x2,x3,x4; y1,y2,y3,y4;z1,z2,z3,z4]);
%     N1=(a1+b1*Dx(i)+c1*Dy(i)+d1*Dz(i))/(6*V);  N2=(a2+b2*Dx(i)+c2*Dy(i)+d2*Dz(i))/(6*V);  N3=(a3+b3*Dx(i)+c3*Dy(i)+d3*Dz(i))/(6*V); N4=(a4+b4*Dx(i)+c4*Dy(i)+d4*Dz(i))/(6*V);
%
%     Mask=ones(1,4);
%     if edgeType(DInd)==1
%         Mask(4)=0;
%     elseif edgeType(DInd)==2
%             Mask(1)=0;
%     elseif edgeType(DInd)==3
%                 Mask(3)=0;
%     else
%                 Mask(2)=0;
%     end
%
%
%
%     Temp=CJ(:,:,i);
%     D1=zeros(Nodes,1);
%     D2=zeros(Nodes,1);
%
%     D1(T(Det.elem(i),:))=Temp(1)*[N1 N2 N3 N4].*Mask;
%     D2(T(Det.elem(i),:))=Temp(2)*[N1 N2 N3 N4].*Mask;
%     D1(D1<1e-12)=0;
%     D2(D2<1e-12)=0;
%     if N==3
%         RxFn(i,:)=[D1',D2'];
%     else
%         D3=zeros(Nodes,1);
%         D3(T(Det.elem(i),:))=Temp(3)*[N1 N2 N3 N4].*Mask;
%         D3(D3<1e-12)=0;
%         if N==5
%             RxFn(i,:)=[D1',D2',D3'];
%         else
%             D4=zeros(Nodes,1);
%             D4(T(Det.elem(i),:))=Temp(4)*[N1 N2 N3 N4].*Mask;
%             D4(D4<1e-12)=0;
%              RxFn(i,:)=[D1',D2',D3',D4'];
%         end
%     end
% end
%
%
% end