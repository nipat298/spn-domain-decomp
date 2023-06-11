% This code works for 3D
function partmat= ddm_interfacing(mesh,NoDomains,varargin)

ElemConnMatrix = [mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3),mesh.tri(:,4)];
NofNodes=size(mesh.nodes,1);
% interface_xyz=['x','x','x'];

t=ElemConnMatrix(:,:);
Nfaces=size(mesh.tri,2);
Nelem=size(mesh.tri,1);

faces = zeros(Nfaces*size(t,1),3); % Forming all possible triangles from the elements 
for ii=1:size(t)
    faces(Nfaces*(ii-1)+1,:) = t(ii,[1,2,3]);
    faces(Nfaces*(ii-1)+2,:) = t(ii,[1,2,4]);
    faces(Nfaces*(ii-1)+3,:) = t(ii,[1,3,4]);
    faces(Nfaces*(ii-1)+4,:) = t(ii,[2,3,4]);
end

% elem_no=zeros(length(ElemConnMatrix(:,1)),Nfaces);
% for i=1:length(ElemConnMatrix(:,1))
%     elem_no(i,:)=[i i i i];
% end
% elem=reshape((elem_no)',(4*384000),1);

elem=(1:Nelem)';
elem=repmat(elem,1,Nfaces);

% for ii=1:length(interface_xyz)
%     switch(interface_xyz(ii))
%         case 'x'

nodes_x12=find(mesh.nodes(:,1)==-0.5000); %Nodes lying on x=-0.5 plane
nodes_x23=find(mesh.nodes(:,1)==0); %Nodes lying on x=0 plane
nodes_x34=find(mesh.nodes(:,1)==0.5000); %Nodes lying on x=0.5 plane


interface_nodes=[nodes_x12 nodes_x23 nodes_x34];


ind1=ismember(faces,nodes_x12); 
% Considering only those triangles from the complete list of triangles
% which have all three nodes lying on x=-0.5 plane
tri_x12_ind= find((ind1(:,1)==1 & ind1(:,2)==1 & ind1(:,3)==1)==1); 
tri_x12=faces(tri_x12_ind,:);
elem_x12=elem(tri_x12_ind);


ind2=ismember(faces,nodes_x23); 
% Considering only those triangles from the complete list of triangles
% which have all three nodes lying on x=0 plane
tri_x23_ind= find((ind2(:,1)==1 & ind2(:,2)==1 & ind2(:,3)==1)==1); 
tri_x23=faces(tri_x23_ind,:);
elem_x23=elem(tri_x23_ind);


ind3=ismember(faces,nodes_x34); 
% Considering only those triangles from the complete list of triangles
% which have all three nodes lying on x=0.5 plane
tri_x34_ind= find((ind3(:,1)==1 & ind3(:,2)==1 & ind3(:,3)==1)==1); 
tri_x34=faces(tri_x34_ind,:);
elem_x34=elem(tri_x34_ind);

n1=[0];
n2=[0];
n3=[0];
n4=[0];

for i=1:NofNodes
    if mesh.nodes(i,1)>=-1 & mesh.nodes(i,1)<-0.5000
        n1=[n1;i];
    elseif mesh.nodes(i,1)>-0.5000 & mesh.nodes(i,1)<0
        n2=[n2;i];
    elseif mesh.nodes(i,1)>0 & mesh.nodes(i,1)<0.5000
        n3=[n3;i];
    elseif mesh.nodes(i,1)>0.5000 & mesh.nodes(i,1)<=1
        n4=[n4;i];
    end
end

n1=unique([n1;nodes_x12]);
n2=unique([n2;nodes_x12;nodes_x23]);
n3=unique([n3;nodes_x23;nodes_x34]);
n4=unique([n4;nodes_x34]);

n1(1,:)=[];
n2(1,:)=[];
n3(1,:)=[];
n4(1,:)=[];

partmat(1).nodes=n1;
partmat(2).nodes=n2;
partmat(3).nodes=n3;
partmat(4).nodes=n4;

partmat(1).interface_g{2}=nodes_x12;
partmat(2).interface_g{1}=nodes_x12;
partmat(2).interface_g{3}=nodes_x23;
partmat(3).interface_g{2}=nodes_x23;
partmat(3).interface_g{4}=nodes_x34;
partmat(4).interface_g{3}=nodes_x34;

for i=1:NoDomains
    for j=1:NoDomains
        if j==(i-1) || j==(i+1)
            continue
        else 
            partmat(i).interface_g{j}=[0];
        end
    end
end


for i=1:NoDomains
    if i==1
        partmat(i).inter_tot=interface_nodes(:,i);
    elseif i==NoDomains
        partmat(i).inter_tot=interface_nodes(:,(i-1));
    else
        partmat(i).inter_tot=unique([interface_nodes(:,(i-1));interface_nodes(:,i)]);
    end
end

% partmat(1).interface_g(2)=nodes_x0;
% partmat(2).interface_g{1}=nodes_x0;


for i=1:NoDomains
    NC_sub=length(partmat(i).nodes);
    partmat(i).SD_nodes=[(1:NC_sub)' partmat(i).nodes];
end

for i=1:NoDomains
    for j=1:NoDomains
%         if j==(i-1) || j==(i+1)
            ismem=ismember(partmat(i).nodes,partmat(i).interface_g{j});
            partmat(i).int_local{j}=find(ismem);
%         else
%              continue
%            
%         end
    end
end

% Forming ElemConnMatrix for subdomains
partmat=ECM_3D(partmat,ElemConnMatrix,NoDomains);
 
%Classifying original boundary and interface triangles for each subdomain
partmat=boundtri_bdelem_3D(partmat,NoDomains);
for i=1:NoDomains
    tt1=ismember(partmat(i).boundtri,partmat(i).inter_tot);
    tt2=find(tt1(:,1)==1 & tt1(:,2)==1 & tt1(:,3)==1);
    
    partmat(i).art_bound=partmat(i).boundtri(tt2,:);
    partmat(i).artb_elem=partmat(i).belem(tt2,:);
    
    partmat(i).boundtri(tt2,:)=[];
    partmat(i).belem(tt2,:)=[];
    
end
