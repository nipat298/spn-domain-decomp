function Svec=GetSource2D_internal(N,P,T,Src)
%Source vector for internal isotropic source embedded in the medium
%Updated on 18 february 2020 - now includes missing area multiplier for
% source element.
% updated on 4 decemeber 2020 - now supports extended sources
Nodes=size(P,1);
D = size(T,2); % Connectivity
% NofSources = length(Src);
NofSrcElems= length(Src.elem);
S=[];

Nmat = ones(NofSrcElems,D);
Nmat = Nmat.*ismember(T(Src.elem,:),Src.node);
Nmat =(Nmat./repmat(sum(Nmat,2),1,D));


for s = 1:NofSrcElems
    S1 = zeros(Nodes,1);
    S1(T(Src.elem(s),:)) = Src.Strength*Nmat(s,:);
    if N==3
        S=[S,[S1;(-2/3)*S1]];
    else if N==5
            S=[S,[S1;-(2/3)*S1; (8/15)*S1]];
        else
            S=[S,[S1;-(2/3)*S1; (8/15)*S1; -(16/35)*S1]];
        end
    end
end

Svec = sum(S,2);

% for ss=1:NofSources
%     
%     if isempty(Src(ss).beam)
%         Sx=Src(ss).Loc(:,1); Sy=Src(ss).Loc(:,2);
%         NofS_elem=size(Src(ss).elem,1);
%         
%         x1=P(T(Src(ss).elem,1),1); x2=P(T(Src(ss).elem,2),1); x3=P(T(Src(ss).elem,3),1);
%         y1=P(T(Src(ss).elem,1),2); y2=P(T(Src(ss).elem,2),2); y3=P(T(Src(ss).elem,3),2);
%         
%         selem = Src(ss).elem;
%         
%         
%     else
%         Sx=Src(ss).beam.loc(:,1); Sy=Src(ss).beam.loc(:,2);
%         NofS_elem=size(Src(ss).beam.elem,1);
%         
%         x1=P(T(Src(ss).beam.elem,1),1); x2=P(T(Src(ss).beam.elem,2),1); x3=P(T(Src(ss).beam.elem,3),1);
%         y1=P(T(Src(ss).beam.elem,1),2); y2=P(T(Src(ss).beam.elem,2),2); y3=P(T(Src(ss).beam.elem,3),2);
%         selem = Src(ss).beam.elem;
%     end
%     
%     a1=x2.*y3-y2.*x3;   a2=x3.*y1-y3.*x1;   a3=x1.*y2-y1.*x2;
%     b1=y2-y3;           b2=y3-y1;           b3=y1-y2;
%     c1=x3-x2;           c2=x1-x3;           c3=x2-x1;
%     A=0.5*(b1.*c2-b2.*c1);
%     
%     N1=(a1+b1.*Sx+c1.*Sy)./(2*A);  N2=(a2+b2.*Sx+c2.*Sy)./(2*A);  N3=(a3+b3.*Sx+c3.*Sy)./(2*A);
%     N1=N1.*A; N2=N2.*A; N3=N3.*A;
%     for i=1:NofS_elem
%         S1=zeros(Nodes,1);
%         S1(T(selem(i),:))=Src.Strength*[N1(i),N2(i),N3(i)];
%         if N==3
%             S=[S,[S1;(-2/3)*S1]];
%         else if N==5
%                 S=[S,[S1;-(2/3)*S1; (8/15)*S1]];
%             else
%                 S=[S,[S1;-(2/3)*S1; (8/15)*S1; -(16/35)*S1]];
%             end
%         end
%     end
%     Svec(:,ss)=sum(S,2);
%     
% end
% trial for distributed source as in Kienle et al.
% dist_2_nodes = sqrt((Src.Loc(1) - P(:,1)).^2 + (Src.Loc(2) - P(:,2)).^2);
% Svec(1:Nodes,1) = 1./(4*pi*dist_2_nodes.^2);
% Svec(Nodes+1:end,1) = (-2/3)*Svec(1:Nodes,1); 
