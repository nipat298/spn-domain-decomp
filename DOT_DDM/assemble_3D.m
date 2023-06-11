function partmat=assemble_3D(DT,partmat,NodeCoord,N,D_tot,muat_tot,Source_loc,mod_freq,speed_in_med,A)

NC=length(NodeCoord(:,1));
nodes_tot=1:length(NodeCoord);

for i=1:N
    
   elem_len=length(partmat(i).ECM(:,1));
    D=D_tot(:,:,partmat(i).ECM(:,5));
    muat=muat_tot(:,:,partmat(i).ECM(:,5));
    
    %% Source Calculation
sourceelem=pointLocation(DT,Source_loc); 
elemindex=find(partmat(i).ECM(:,5)==sourceelem);

[ID,d]=nearestNeighbor(DT,Source_loc);
s=length(ID); %s is the number of sources

source=zeros(length(NodeCoord),1);
 
if isempty(elemindex)==0
     n1=partmat(i).ECM(elemindex,1); 
     n2=partmat(i).ECM(elemindex,2);
     n3=partmat(i).ECM(elemindex,3);
     n4=partmat(i).ECM(elemindex,4);

      x1s=NodeCoord(n1,1); x2s=NodeCoord(n2,1); x3s=NodeCoord(n3,1); x4s=NodeCoord(n4,1);
      y1s=NodeCoord(n1,2); y2s=NodeCoord(n2,2); y3s=NodeCoord(n3,2); y4s=NodeCoord(n4,2);
      z1s=NodeCoord(n1,3); z2s=NodeCoord(n2,3); z3s=NodeCoord(n3,3); z4s=NodeCoord(n4,3);

    %  elemvol_source=(1/6)*((det([x2s x3s x4s; y2s y3s y4s; z2s z3s z4s]))-(det([x1s x3s x4s; y1s y3s y4s; z1s z3s z4s]))+(det([x1s x2s x4s; y1s y2s y4s; z1s z2s z4s]))-(det([x1s x2s x3s; y1s y2s y3s; z1s z2s z3s])));
%       det1=det([x2s x3s x4s; y2s y3s y4s; z2s z3s z4s]);
%       det2=det([x1s x3s x4s; y1s y3s y4s; z1s z3s z4s]);
%       det3=det([x1s x2s x4s; y1s y2s y4s; z1s z2s z4s]);
%       det4=det([x1s x2s x3s; y1s y2s y3s; z1s z2s z3s]);
%       elemvol_source=(1/6)*(det1-det2+det3-det4);

        a1s=x2s.*(y3s.*z4s-y4s.*z3s) - x3s.*(y2s.*z4s-y4s.*z2s) + x4s.*(y2s.*z3s-y3s.*z2s);
        a2s=-(x1s.*(y3s.*z4s-y4s.*z3s) - x3s.*(y1s.*z4s-y4s.*z1s) + x4s.*(y1s.*z3s-y3s.*z1s));
        a3s=x1s.*(y2s.*z4s-y4s.*z2s) - x2s.*(y1s.*z4s-y4s.*z1s) + x4s.*(y1s.*z2s-y2s.*z1s);
        a4s=-(x1s.*(y2s.*z3s-y3s.*z2s) - x2s.*(y1s.*z3s-y3s.*z1s) + x3s.*(y1s.*z2s-y2s.*z1s));

        b1s=-((y3s.*z4s-y4s.*z3s) - (y2s.*z4s-y4s.*z2s) + (y2s.*z3s-y3s.*z2s));
        b2s=(y3s.*z4s-y4s.*z3s) - (y1s.*z4s-y4s.*z1s) + (y1s.*z3s-y3s.*z1s);
        b3s=-((y2s.*z4s-y4s.*z2s) - (y1s.*z4s-y4s.*z1s) + (y1s.*z2s-y2s.*z1s));
        b4s=(y2s.*z3s-y3s.*z2s) - (y1s.*z3s-y3s.*z1s) + (y1s.*z2s-y2s.*z1s);

        c1s=(x3s.*z4s-x4s.*z3s) - (x2s.*z4s-x4s.*z2s) + (x2s.*z3s-x3s.*z2s);
        c2s=-((x3s.*z4s-x4s.*z3s) - (x1s.*z4s-x4s.*z1s) + (x1s.*z3s-x3s.*z1s));
        c3s=(x2s.*z4s-x4s.*z2s) - (x1s.*z4s-x4s.*z1s) + (x1s.*z2s-x2s.*z1s);
        c4s=-((x2s.*z3s-x3s.*z2s) - (x1s.*z3s-x3s.*z1s) + (x1s.*z2s-x2s.*z1s));

        d1s=(y3s.*x4s-y4s.*x3s) - (y2s.*x4s-y4s.*x2s) + (y2s.*x3s-y3s.*x2s);
        d2s=-((y3s.*x4s-y4s.*x3s) - (y1s.*x4s-y4s.*x1s) + (y1s.*x3s-y3s.*x1s));
        d3s=(y2s.*x4s-y4s.*x2s) - (y1s.*x4s-y4s.*x1s) + (y1s.*x2s-y2s.*x1s);
        d4s=-((y2s.*x3s-y3s.*x2s) - (y1s.*x3s-y3s.*x1s) + (y1s.*x2s-y2s.*x1s));

        vol_elem=(1/6)*det([1,1,1,1; x1s,x2s,x3s,x4s; y1s,y2s,y3s,y4s; z1s,z2s,z3s,z4s]);

      xcoord=Source_loc(1,1);
      ycoord=Source_loc(1,2);
      zcoord=Source_loc(1,3);
      
      source(n1,1)=(a1s+ b1s.*xcoord+ c1s.*ycoord+ d1s.*zcoord)./(6*vol_elem);
      source(n2,1)=(a2s+ b2s.*xcoord+ c2s.*ycoord+ d2s.*zcoord)./(6*vol_elem);
      source(n3,1)=(a3s +b3s.*xcoord+ c3s.*ycoord+ d3s.*zcoord)./(6*vol_elem);
      source(n4,1)=(a4s +b4s.*xcoord+ c4s.*ycoord+ d4s.*zcoord)./(6*vol_elem);


%       source(n1,1)=((det([x2s x3s x4s; y2s y3s y4s; z2s z3s z4s]))-(det([xcoord x3s x4s; ycoord y3s y4s; zcoord z3s z4s]))+(det([xcoord x2s x4s; ycoord y2s y4s; zcoord z2s z4s]))-(det([xcoord x2s x3s; ycoord y2s y3s; zcoord z2s z3s])))/abs(elemvol_source);
%       source(n2,1)=((det([xcoord x3s x4s; ycoord y3s y4s; zcoord z3s z4s]))-(det([x1s x3s x4s; y1s y3s y4s; z1s z3s z4s]))+(det([x1s xcoord x4s; y1s ycoord y4s; z1s zcoord z4s]))-(det([x1s xcoord x3s; y1s ycoord y3s; z1s zcoord z3s])))/abs(elemvol_source);
%       source(n3,1)=((det([x2s xcoord x4s; y2s ycoord y4s; z2s zcoord z4s]))-(det([x1s xcoord x4s; y1s ycoord y4s; z1s zcoord z4s]))+(det([x1s x2s x4s; y1s y2s y4s; z1s z2s z4s]))-(det([x1s x2s xcoord; y1s y2s ycoord; z1s z2s zcoord])))/abs(elemvol_source);
%       source(n4,1)=((det([x2s x3s xcoord; y2s y3s ycoord; z2s z3s zcoord]))-(det([x1s x3s xcoord; y1s y3s ycoord; z1s z3s zcoord]))+(det([x1s x2s xcoord; y1s y2s ycoord; z1s z2s zcoord]))-(det([x1s x2s x3s; y1s y2s y3s; z1s z2s z3s])))/abs(elemvol_source);

end
      partmat(i).source=source;


%% Matrix Assembly

    a1=zeros(1,1,elem_len);
    a2=zeros(1,1,elem_len);
    a3=zeros(1,1,elem_len);
    a4=zeros(1,1,elem_len);
    
    b1=zeros(1,1,elem_len);
    b2=zeros(1,1,elem_len);
    b3=zeros(1,1,elem_len);
    b4=zeros(1,1,elem_len);
    
    c1=zeros(1,1,elem_len);
    c2=zeros(1,1,elem_len);
    c3=zeros(1,1,elem_len);
    c4=zeros(1,1,elem_len);
    
    d1=zeros(1,1,elem_len);
    d2=zeros(1,1,elem_len);
    d3=zeros(1,1,elem_len);
    d4=zeros(1,1,elem_len);
    
    part1=zeros(1,1,elem_len);
    part2=zeros(1,1,elem_len);
    part3=zeros(1,1,elem_len);
    elem_vol=zeros(1,1,elem_len);



partmat(i).ECM(:,5)=[];


x1=NodeCoord(partmat(i).ECM(:,1),1);
x2=NodeCoord(partmat(i).ECM(:,2),1);
x3=NodeCoord(partmat(i).ECM(:,3),1);
x4=NodeCoord(partmat(i).ECM(:,4),1);


y1=NodeCoord(partmat(i).ECM(:,1),2);
y2=NodeCoord(partmat(i).ECM(:,2),2);
y3=NodeCoord(partmat(i).ECM(:,3),2);
y4=NodeCoord(partmat(i).ECM(:,4),2);

z1=NodeCoord(partmat(i).ECM(:,1),3);
z2=NodeCoord(partmat(i).ECM(:,2),3);
z3=NodeCoord(partmat(i).ECM(:,3),3);
z4=NodeCoord(partmat(i).ECM(:,4),3);

a1(1,1,:) = x2.*(y3.*z4 - y4.*z3) + x3.*(y4.*z2 - y2.*z4) + x4.*(y2.*z3 - y3.*z2);
a2(1,1,:) = x1.*(y4.*z3 - y3.*z4) + x3.*(y1.*z4 - y4.*z1) + x4.*(y3.*z1 - y1.*z3);
a3(1,1,:) = x1.*(y2.*z4 - y4.*z2) + x2.*(y4.*z1 - y1.*z4) + x4.*(y1.*z2 - y2.*z1);
a4(1,1,:) = x1.*(y3.*z2 - y2.*z3) + x2.*(y1.*z3 - y3.*z1) + x3.*(y2.*z1 - y1.*z2);

b1(1,1,:) = y2.*(z4-z3) + y3.*(z2-z4) + y4.*(z3-z2);
b2(1,1,:) = y1.*(z3-z4) + y3.*(z4-z1) + y4.*(z1-z3);
b3(1,1,:) = y1.*(z4-z2) + y2.*(z1-z4) + y4.*(z2-z1);
b4(1,1,:) = y1.*(z2-z3) + y2.*(z3-z1) + y3.*(z1-z2);

c1(1,1,:) = z2.*(x4-x3) + z3.*(x2-x4) + z4.*(x3-x2);
c2(1,1,:) = z1.*(x3-x4) + z3.*(x4-x1) + z4.*(x1-x3);
c3(1,1,:) = z1.*(x4-x2) + z2.*(x1-x4) + z4.*(x2-x1);
c4(1,1,:) = z1.*(x2-x3) + z2.*(x3-x1) + z3.*(x1-x2);

d1(1,1,:) = x2.*(y4-y3) + x3.*(y2-y4) + x4.*(y3-y2);
d2(1,1,:) = x1.*(y3-y4) + x3.*(y4-y1) + x4.*(y1-y3);
d3(1,1,:) = x1.*(y4-y2) + x2.*(y1-y4) + x4.*(y2-y1);
d4(1,1,:) = x1.*(y2-y3) + x2.*(y3-y1) + x3.*(y1-y2);

%%% Determinant of Jacobian

part1(1,1,:)=(x2-x1).*((y3-y1).*(z4-z1)-(z3-z1).*(y4-y1));
part2(1,1,:)=(x3-x1).*((y4-y1).*(z2-z1)-(y2-y1).*(z4-z1));
part3(1,1,:)=(x4-x1).*((y2-y1).*(z3-z1)-(y3-y1).*(z2-z1));
elem_vol(1,1,:)=(1/6)*(part1+part2+part3);


%%% Stiffness matrix
Ks=zeros(4,4,elem_len);

% % Ks(1,1,:)=(b1.^2+ c1.^2+ d1.^2)./(6*elem_vol);
% % Ks(1,2,:)=(b1.*b2+ c1.*c2+ d1.*d2)./(6*elem_vol);
% % Ks(1,3,:)=(b1.*b3+ c1.*c3+ d1.*d3)./(6*elem_vol);
% % Ks(1,4,:)=(b1.*b4+ c1.*c4+ d1.*d4)./(6*elem_vol);
% % 
% % Ks(2,1,:)=Ks(1,2,:);
% % Ks(2,2,:)=(b2.^2+ c2.^2+ d2.^2)./(6*elem_vol);
% % Ks(2,3,:)=(b2.*b3+ c2.*c3+ d2.*d3)./(6*elem_vol);
% % Ks(2,4,:)=(b2.*b4+ c2.*c4+ d2.*d4)./(6*elem_vol);
% % 
% % Ks(3,1,:)=Ks(1,3,:);
% % Ks(3,2,:)=Ks(2,3,:);
% % Ks(3,3,:)=(b3.^2+ c3.^2+ d3.^2)./(6*elem_vol);
% % Ks(3,4,:)=(b3.*b4+ c3.*c4+ d3.*d4)./(6*elem_vol);
% % 
% % Ks(4,1,:)=Ks(1,4,:);
% % Ks(4,2,:)=Ks(2,4,:);
% % Ks(4,3,:)=Ks(3,4,:);
% % Ks(4,4,:)=(b4.^2+ c4.^2+ d4.^2)./(6*elem_vol);



Ks(1,1,:)=((b1(1,1,:).^2)+(c1(1,1,:).^2)+(d1(1,1,:).^2))./(36*elem_vol);
Ks(1,2,:)=((b1(1,1,:).*b2(1,1,:))+(c1(1,1,:).*c2(1,1,:))+(d1(1,1,:).*d2(1,1,:)))./(36*elem_vol);
Ks(1,3,:)=((b1(1,1,:).*b3(1,1,:))+(c1(1,1,:).*c3(1,1,:))+(d1(1,1,:).*d3(1,1,:)))./(36*elem_vol);
Ks(1,4,:)=((b1(1,1,:).*b4(1,1,:))+(c1(1,1,:).*c4(1,1,:))+(d1(1,1,:).*d4(1,1,:)))./(36*elem_vol);

Ks(2,1,:)=Ks(1,2,:);
Ks(2,2,:)=((b2(1,1,:).^2)+(c2(1,1,:).^2)+(d2(1,1,:).^2))./(36*elem_vol);
Ks(2,3,:)=((b2(1,1,:).*b3(1,1,:))+(c2(1,1,:).*c3(1,1,:))+(d2(1,1,:).*d3(1,1,:)))./(36*elem_vol);
Ks(2,4,:)=((b2(1,1,:).*b4(1,1,:))+(c2(1,1,:).*c4(1,1,:))+(d2(1,1,:).*d4(1,1,:)))./(36*elem_vol);

Ks(3,1,:)=Ks(1,3,:);
Ks(3,2,:)=Ks(2,3,:);
Ks(3,3,:)=((b3(1,1,:).^2)+(c3(1,1,:).^2)+(d3(1,1,:).^2))./(36*elem_vol);
Ks(3,4,:)=((b3(1,1,:).*b4(1,1,:))+(c3(1,1,:).*c4(1,1,:))+(d3(1,1,:).*d4(1,1,:)))./(36*elem_vol);

Ks(4,1,:)=Ks(1,4,:);
Ks(4,2,:)=Ks(2,4,:);
Ks(4,3,:)=Ks(3,4,:);
Ks(4,4,:)=((b4(1,1,:).^2)+(c4(1,1,:).^2)+(d4(1,1,:).^2))./(36*elem_vol);
Kstiff=repmat(D(1,1,:),4,4).*Ks;


%%% Mass matrix
Km=zeros(4,4,elem_len);
Km=[elem_vol/10, elem_vol/20, elem_vol/20, elem_vol/20
    elem_vol/20, elem_vol/10, elem_vol/20, elem_vol/20
    elem_vol/20, elem_vol/20, elem_vol/10, elem_vol/20
    elem_vol/20, elem_vol/20, elem_vol/20, elem_vol/10];
Kmass=repmat(muat(1,1,:),4,4).*Km+ ((1i*mod_freq)/speed_in_med).*Km;



%%% Boundary matrix
tri_area=zeros(1,1,elem_len); 
% In the following eqn, we calculate the area of each boundary triangle and
% assign it to the element corresponding to it
    
bound_tri(i).a=sqrt((NodeCoord(partmat(i).boundtri(:,1),1)-NodeCoord(partmat(i).boundtri(:,2),1)).^2+(NodeCoord(partmat(i).boundtri(:,1),2)-NodeCoord(partmat(i).boundtri(:,2),2)).^2+(NodeCoord(partmat(i).boundtri(:,1),3)-NodeCoord(partmat(i).boundtri(:,2),3)).^2);
bound_tri(i).b=sqrt((NodeCoord(partmat(i).boundtri(:,2),1)-NodeCoord(partmat(i).boundtri(:,3),1)).^2+(NodeCoord(partmat(i).boundtri(:,2),2)-NodeCoord(partmat(i).boundtri(:,3),2)).^2+(NodeCoord(partmat(i).boundtri(:,2),3)-NodeCoord(partmat(i).boundtri(:,3),3)).^2);
bound_tri(i).c=sqrt((NodeCoord(partmat(i).boundtri(:,1),1)-NodeCoord(partmat(i).boundtri(:,3),1)).^2+(NodeCoord(partmat(i).boundtri(:,1),2)-NodeCoord(partmat(i).boundtri(:,3),2)).^2+(NodeCoord(partmat(i).boundtri(:,1),3)-NodeCoord(partmat(i).boundtri(:,3),3)).^2);
bound_tri(i).s=(bound_tri(i).a+bound_tri(i).b+bound_tri(i).c)/2;

tri_area(1,1,partmat(i).belem(:,1))=sqrt((bound_tri(i).s).*(bound_tri(i).s-bound_tri(i).a).*(bound_tri(i).s-bound_tri(i).b).*(bound_tri(i).s-bound_tri(i).c));

f_123=[1/6 1/12 1/12 0; 1/12 1/6 1/12 0; 1/12 1/12 1/6 0; 0 0 0 0];
f_124=[1/6 1/12 0 1/12; 1/12 1/6 0 1/12; 0 0 0 0; 1/12 1/12 0 1/6];
f_134=[1/6 0 1/12 1/12; 0 0 0 0; 1/12 0 1/6 1/12; 1/12 0 1/12 1/6];
f_234=[0 0 0 0; 0 1/6 1/12 1/12; 0 1/12 1/6 1/12; 0 1/12 1/12 1/6];


Kb=zeros(4,4,elem_len);

i1=find(partmat(i).belem(:,2)==1);
Kb(:,:,partmat(i).belem(i1,1))=repmat(f_123,[1 1 size(i1,1)]).*repmat(tri_area(partmat(i).belem(i1,1)),[4,4]);
i1=find(partmat(i).belem(:,2)==2);
Kb(:,:,partmat(i).belem(i1,1))=repmat(f_234,[1 1 size(i1,1)]).*repmat(tri_area(partmat(i).belem(i1,1)),[4,4]);
i1=find(partmat(i).belem(:,2)==3);
Kb(:,:,partmat(i).belem(i1,1))=repmat(f_124,[1 1 size(i1,1)]).*repmat(tri_area(partmat(i).belem(i1,1)),[4,4]);
i1=find(partmat(i).belem(:,2)==4);
Kb(:,:,partmat(i).belem(i1,1))=repmat(f_134,[1 1 size(i1,1)]).*repmat(tri_area(partmat(i).belem(i1,1)),[4,4]);
Kbound=(0.5/A).*Kb;


Ti = reshape((partmat(i).ECM)',[4,1,elem_len]);
Tj = reshape((partmat(i).ECM)',[1,4,elem_len]);
I = repmat(Ti,1,4);
J = repmat(Tj,4,1);

Kstiff_final = sparse(I(:),J(:),Kstiff(:),NC,NC);
Kmass_final = sparse(I(:),J(:),Kmass(:),NC,NC);
Kbound_final = sparse(I(:),J(:),Kbound(:),NC,NC);

% Assembling of global sparse matrix
G=Kstiff_final+Kmass_final+Kbound_final;

partmat(i).Gsub=G;


sub_ind=setdiff((nodes_tot)',partmat(i).nodes);
% int_SD_index=find(ismember(inter_tri,partmat(i).nodes)==1);
% inter_SD=inter_tri(int_SD_index);

% partmat(i).Gsub(inter_SD,:)=0;
% partmat(i).Gsub(:,inter_SD)=0;

partmat(i).Gsub(sub_ind,:)=[];
partmat(i).Gsub(:,sub_ind)=[];
partmat(i).source(sub_ind,:)=[];





end







