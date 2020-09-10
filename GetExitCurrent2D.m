function RxFn=GetExitCurrent2D(N,CJ,P,T,edgeElem,edgeType,Det)
%Compute the CJ term in equation (2)
d=numel(Det.elem);
Dx=Det.Loc(:,1);Dy=Det.Loc(:,2);
%     T=[DT(:,1),DT(:,2),DT(:,3)];
Nodes=size(P,1);
%     Elem=size(T,1);

%     [~,~,e_n,~]=GetEdgeLength(P,T);
for i=1:d
    DInd=edgeElem(:,1)==Det.elem(i);  %Index
    
    x1=P(T(Det.elem(i),1),1); x2=P(T(Det.elem(i),2),1); x3=P(T(Det.elem(i),3),1);
    y1=P(T(Det.elem(i),1),2); y2=P(T(Det.elem(i),2),2); y3=P(T(Det.elem(i),3),2);
    
    a1=x2.*y3-y2.*x3;   a2=x3.*y1-y3.*x1;   a3=x1.*y2-y1.*x2;
    b1=y2-y3;           b2=y3-y1;           b3=y1-y2;
    c1=x3-x2;           c2=x1-x3;           c3=x2-x1;
    A=0.5*(b1.*c2-b2.*c1);
    
    N1=(a1+b1*Dx(i)+c1*Dy(i))/(2*A);  N2=(a2+b2*Dx(i)+c2*Dy(i))/(2*A);  N3=(a3+b3*Dx(i)+c3*Dy(i))/(2*A);
    
    Mask=ones(1,3);
    if edgeType(DInd)==1
        Mask(3)=0;
%        Mask(1:2)=sqrt((x1-x2).^2 + (y1-y2).^2);
    else if edgeType(DInd)==2
            Mask(1)=0;
%        Mask(2:3)=sqrt((x3-x2).^2 + (y3-y2).^2);
        else
            Mask(2)=0;
%            Mask([1,3])=sqrt((x1-x3).^2 + (y1-y3).^2);
        end
    end
    
    Temp=CJ(:,:,i);
    D1=zeros(Nodes,1);
    D1(T(Det.elem(i),:))=Temp(1)*[N1 N2 N3].*Mask;
    D1(D1<1e-12)=0;
    if N==1
        RxFn(i,:)=D1';
    else 
         D2=zeros(Nodes,1);
         D2(T(Det.elem(i),:))=Temp(2)*[N1 N2 N3].*Mask;
         D2(D2<1e-12)=0;
         if N==3
            RxFn(i,:)=[D1',D2'];
         else
            D3=zeros(Nodes,1);
            D3(T(Det.elem(i),:))=Temp(3)*[N1 N2 N3].*Mask;
            if N==5
                RxFn(i,:)=[D1',D2',D3'];
            else
                D4=zeros(Nodes,1);
                D4(T(Det.elem(i),:))=Temp(4)*[N1 N2 N3].*Mask;
                RxFn(i,:)=[D1',D2',D3',D4'];
            end
        end
    end
end

end