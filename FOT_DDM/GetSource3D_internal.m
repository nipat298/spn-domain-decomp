function Svec=GetSource3D_internal(N,P,T,Src)
% Modified on 5th December 2020 - supports extended sources
%Source vector for internal isotropic source embedded in the medium
%Updated on 18 february 2020 - now includes missing area multiplier for

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

% source element.
% Nodes=size(P,1);
% NofSources = length(Src);
% S=[];
% Svec = zeros((N+1)/2*Nodes,NofSources);
% 
% for ss=1:NofSources
%     
%     if isempty(Src(ss).beam)
%         Sx=Src(ss).Loc(:,1); Sy=Src(ss).Loc(:,2); Sz=Src(ss).Loc(:,3);
%         NofS_elem=size(Src(ss).elem,1);
%         
%         x1=P(T(Src(ss).elem,1),1); x2=P(T(Src(ss).elem,2),1); x3=P(T(Src(ss).elem,3),1); x4=P(T(Src(ss).elem,4),1);
%         y1=P(T(Src(ss).elem,1),2); y2=P(T(Src(ss).elem,2),2); y3=P(T(Src(ss).elem,3),2); y4=P(T(Src(ss).elem,4),2);
%         z1=P(T(Src(ss).elem,1),3); z2=P(T(Src(ss).elem,2),3); z3=P(T(Src(ss).elem,3),3); z4=P(T(Src(ss).elem,4),3);
%         
%         selem = Src(ss).elem;
%         
%         
%     else
%         Sx=Src(ss).beam.loc(:,1); Sy=Src(ss).beam.loc(:,2); Sz=Src(ss).beam.loc(:,3);
%         NofS_elem=size(Src(ss).beam.elem,1);
%         
%         x1=P(T(Src(ss).beam.elem,1),1); x2=P(T(Src(ss).beam.elem,2),1); x3=P(T(Src(ss).beam.elem,3),1); x4=P(T(Src(ss).beam.elem,4),1);
%         y1=P(T(Src(ss).beam.elem,1),2); y2=P(T(Src(ss).beam.elem,2),2); y3=P(T(Src(ss).beam.elem,3),2); y4=P(T(Src(ss).beam.elem,4),2);
%         z1=P(T(Src(ss).beam.elem,1),3); z2=P(T(Src(ss).beam.elem,2),3); z3=P(T(Src(ss).beam.elem,3),3); z4=P(T(Src(ss).beam.elem,4),3);
%         selem = Src(ss).beam.elem;
%     end
%     
%    a1=x2.*(y3.*z4-y4.*z3) - x3.*(y2.*z4-y4.*z2) + x4.*(y2.*z3-y3.*z2);
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
% %  N1=N1.*A; N2=N2.*A; N3=N3.*A;
%     for i=1:NofS_elem
%         S1=zeros(Nodes,1);
%         S1(T(selem(i),:))=Src.Strength*[N1(i),N2(i),N3(i),N4(i)];
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
