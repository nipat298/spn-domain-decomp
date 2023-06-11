%%% Testing the domain decomposition code with manual partitioning without the help of METIS
%%% 
% % % % % % % % % % % % tic
% % % % % % % % % % % % mesh_spacing = 0.05;
% % % % % % % % % % % % xmin = -1; xmax =1; ymin = -1; ymax = 1; zmin = 0; zmax = 2;
% % % % % % % % % % % % [X,Y,Z] = meshgrid((xmin:mesh_spacing:xmax),(ymin:mesh_spacing:ymax),(zmin:mesh_spacing:zmax)); % generates an equispaced grid of nodes
% % % % % % % % % % % % mesh.nodes = [X(:),Y(:),Z(:)];     % x,y coordinates of the nodes
% % % % % % % % % % % % DT =  delaunayTriangulation(mesh.nodes);
% % % % % % % % % % % % mesh.tri = [DT(:,1),DT(:,2),DT(:,3),DT(:,4)];
% % % % % % % % % % % % NC=length(mesh.nodes);
% % % % % % % % % % % % N=4;
% % % % % % % % % % % % 
% % % % % % % % % % % % t=mesh.tri(:,:);
% % % % % % % % % % % % faces = zeros(4*size(t,1),3); % Forming all possible triangles from the elements 
% % % % % % % % % % % % for ii=1:size(t)
% % % % % % % % % % % %     faces(4*(ii-1)+1,:) = t(ii,[1,2,3]);
% % % % % % % % % % % %     faces(4*(ii-1)+2,:) = t(ii,[1,2,4]);
% % % % % % % % % % % %     faces(4*(ii-1)+3,:) = t(ii,[1,3,4]);
% % % % % % % % % % % %     faces(4*(ii-1)+4,:) = t(ii,[2,3,4]);
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % elem_no=zeros(length(mesh.tri(:,1)),4);
% % % % % % % % % % % % for i=1:length(mesh.tri(:,1))
% % % % % % % % % % % %     elem_no(i,:)=[i i i i];
% % % % % % % % % % % % end
% % % % % % % % % % % % elem=reshape((elem_no)',(4*384000),1);
% % % % % % % % % % % % 
% % % % % % % % % % % % nodes_x12=find(mesh.nodes(:,1)==-0.5000); %Nodes lying on x=-0.5 plane
% % % % % % % % % % % % nodes_x23=find(mesh.nodes(:,1)==0); %Nodes lying on x=0 plane
% % % % % % % % % % % % nodes_x34=find(mesh.nodes(:,1)==0.5000); %Nodes lying on x=0.5 plane
% % % % % % % % % % % % 
% % % % % % % % % % % % interface_nodes=[nodes_x12 nodes_x23 nodes_x34];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % % % % ind1=ismember(faces,nodes_x12); 
% % % % % % % % % % % % % % % % % Considering only those triangles from the complete list of triangles
% % % % % % % % % % % % % % % % % which have all three nodes lying on x=-0.5 plane
% % % % % % % % % % % % % % % % tri_x12_ind= find((ind1(:,1)==1 & ind1(:,2)==1 & ind1(:,3)==1)==1); 
% % % % % % % % % % % % % % % % tri_x12=faces(tri_x12_ind,:);
% % % % % % % % % % % % % % % % elem_x12=elem(tri_x12_ind,:);
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % ind2=ismember(faces,nodes_x23); 
% % % % % % % % % % % % % % % % % Considering only those triangles from the complete list of triangles
% % % % % % % % % % % % % % % % % which have all three nodes lying on x=-0.5 plane
% % % % % % % % % % % % % % % % tri_x23_ind= find((ind2(:,1)==1 & ind2(:,2)==1 & ind2(:,3)==1)==1); 
% % % % % % % % % % % % % % % % tri_x23=faces(tri_x23_ind,:);
% % % % % % % % % % % % % % % % elem_x23=elem(tri_x23_ind,:);
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % ind3=ismember(faces,nodes_x34); 
% % % % % % % % % % % % % % % % % Considering only those triangles from the complete list of triangles
% % % % % % % % % % % % % % % % % which have all three nodes lying on x=-0.5 plane
% % % % % % % % % % % % % % % % tri_x34_ind= find((ind3(:,1)==1 & ind3(:,2)==1 & ind3(:,3)==1)==1); 
% % % % % % % % % % % % % % % % tri_x34=faces(tri_x34_ind,:);
% % % % % % % % % % % % % % % % elem_x34=elem(tri_x34_ind,:);
% % % % % % % % % % % % 
% % % % % % % % % % % % n1=[0];
% % % % % % % % % % % % n2=[0];
% % % % % % % % % % % % n3=[0];
% % % % % % % % % % % % n4=[0];
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:NC
% % % % % % % % % % % %     if mesh.nodes(i,1)>=-1 & mesh.nodes(i,1)<-0.5000
% % % % % % % % % % % %         n1=[n1;i];
% % % % % % % % % % % %     elseif mesh.nodes(i,1)>-0.5000 & mesh.nodes(i,1)<0
% % % % % % % % % % % %         n2=[n2;i];
% % % % % % % % % % % %     elseif mesh.nodes(i,1)>0 & mesh.nodes(i,1)<0.5000
% % % % % % % % % % % %         n3=[n3;i];
% % % % % % % % % % % %     elseif mesh.nodes(i,1)>0.5000 & mesh.nodes(i,1)<=1
% % % % % % % % % % % %         n4=[n4;i];
% % % % % % % % % % % %     end
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % n1=unique([n1;nodes_x12]);
% % % % % % % % % % % % n2=unique([n2;nodes_x12;nodes_x23]);
% % % % % % % % % % % % n3=unique([n3;nodes_x23;nodes_x34]);
% % % % % % % % % % % % n4=unique([n4;nodes_x34]);
% % % % % % % % % % % % 
% % % % % % % % % % % % n1(1,:)=[];
% % % % % % % % % % % % n2(1,:)=[];
% % % % % % % % % % % % n3(1,:)=[];
% % % % % % % % % % % % n4(1,:)=[];
% % % % % % % % % % % % 
% % % % % % % % % % % % partmat(1).nodes=n1;
% % % % % % % % % % % % partmat(2).nodes=n2;
% % % % % % % % % % % % partmat(3).nodes=n3;
% % % % % % % % % % % % partmat(4).nodes=n4;
% % % % % % % % % % % % 
% % % % % % % % % % % % % for i=1:N
% % % % % % % % % % % % %     for j=1:N
% % % % % % % % % % % % %         if j==i
% % % % % % % % % % % % %             partmat(i).interface_g{j}=[];
% % % % % % % % % % % % %         else
% % % % % % % % % % % % %             partmat(i).interface_g{j}=nodes_x0;
% % % % % % % % % % % % %         end
% % % % % % % % % % % % %     end
% % % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % partmat(1).interface_g{2}=nodes_x12;
% % % % % % % % % % % % partmat(2).interface_g{1}=nodes_x12;
% % % % % % % % % % % % partmat(2).interface_g{3}=nodes_x23;
% % % % % % % % % % % % partmat(3).interface_g{2}=nodes_x23;
% % % % % % % % % % % % partmat(3).interface_g{4}=nodes_x34;
% % % % % % % % % % % % partmat(4).interface_g{3}=nodes_x34;
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:N
% % % % % % % % % % % %     for j=1:N
% % % % % % % % % % % %         if j==(i-1) || j==(i+1)
% % % % % % % % % % % %             continue
% % % % % % % % % % % %         else 
% % % % % % % % % % % %             partmat(i).interface_g{j}=[0];
% % % % % % % % % % % %         end
% % % % % % % % % % % %     end
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:N
% % % % % % % % % % % %     if i==1
% % % % % % % % % % % %         partmat(i).inter_tot=interface_nodes(:,i);
% % % % % % % % % % % %     elseif i==N
% % % % % % % % % % % %         partmat(i).inter_tot=interface_nodes(:,(i-1));
% % % % % % % % % % % %     else
% % % % % % % % % % % %         partmat(i).inter_tot=unique([interface_nodes(:,(i-1));interface_nodes(:,i)]);
% % % % % % % % % % % %     end
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % % partmat(1).interface_g(2)=nodes_x0;
% % % % % % % % % % % % % partmat(2).interface_g{1}=nodes_x0;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:N
% % % % % % % % % % % %     NC_sub=length(partmat(i).nodes);
% % % % % % % % % % % %     partmat(i).SD_nodes=[(1:NC_sub)' partmat(i).nodes];
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:N
% % % % % % % % % % % %     for j=1:N
% % % % % % % % % % % % %         if j==(i-1) || j==(i+1)
% % % % % % % % % % % %             ismem=ismember(partmat(i).nodes,partmat(i).interface_g{j});
% % % % % % % % % % % %             partmat(i).int_local{j}=find(ismem);
% % % % % % % % % % % % %         else
% % % % % % % % % % % % %              continue
% % % % % % % % % % % % %            
% % % % % % % % % % % % %         end
% % % % % % % % % % % %     end
% % % % % % % % % % % % end

function partmat=noMETIS_4sub(partmat,N,mesh,alpha,beta,DT,D,muat,Source_loc,mod_freq,speed_in_med,A)


%  [e,e_n]=boundedges_element_3D(mesh.nodes,mesh.tri);
%  ind1=ismember(e,partmat(1).nodes);
%  ind_n1=find(ind1(:,1)==1 & ind1(:,2)==1 & ind1(:,3)==1);
%  partmat(1).boundtri=e(ind_n1,:);
%  partmat(1).belem=e_n(ind_n1,:);
%  
%  ind2=ismember(e,partmat(2).nodes);
%  ind_n2=find(ind2(:,1)==1 & ind2(:,2)==1 & ind2(:,3)==1);
%  partmat(2).boundtri=e(ind_n2,:);
%  partmat(2).belem=e_n(ind_n2,:);
 
% partmat=boundtri_bdelem_3D(partmat,N);
% tt1=ismember(partmat(1).boundtri,nodes_x0);
% tt2=find(tt1(:,1)==1 & tt1(:,2)==1 & tt1(:,3)==1);
% 
% partmat(1).art_bound=partmat(1).boundtri(tt2,:);
% partmat(1).artb_elem=partmat(1).belem(tt2,:);
% 
% partmat(1).boundtri(tt2,:)=[];
% partmat(1).belem(tt2,:)=[];
% 
% tt3=ismember(partmat(2).boundtri,nodes_x0);
% tt4=find(tt3(:,1)==1 & tt3(:,2)==1 & tt3(:,3)==1);

% partmat(2).art_bound=partmat(2).boundtri(tt4,:);
% partmat(2).artb_elem=partmat(2).belem(tt4,:);
% 
% partmat(2).boundtri(tt4,:)=[];
% partmat(2).belem(tt4,:)=[];


 partmat=ECM_3D(partmat,mesh.tri,N);
 
partmat=boundtri_bdelem_3D(partmat,N);
for i=1:N
    tt1=ismember(partmat(i).boundtri,partmat(i).inter_tot);
    tt2=find(tt1(:,1)==1 & tt1(:,2)==1 & tt1(:,3)==1);
    
    partmat(i).art_bound=partmat(i).boundtri(tt2,:);
    partmat(i).artb_elem=partmat(i).belem(tt2,:);
    
    partmat(i).boundtri(tt2,:)=[];
    partmat(i).belem(tt2,:)=[];
end


%  partmat(1).art_bound=tri_x0;
%  partmat(2).art_bound=tri_x0;
%  
%  partmat(1).artb_elem=elem_x0;
%  partmat(2).artb_elem=elem_x0;
 
  partmat=assemble_3D(DT,partmat,mesh.nodes,N,D,muat,Source_loc,mod_freq,speed_in_med,A);
  
%    beta=1/(mesh_spacing^2);
%  alpha=1/(mesh_spacing);
% % 
% %    beta=1/((3*mesh_spacing)^2);
% %   alpha=1/(0.1*mesh_spacing);
 partmat=interface_calc_3D(partmat,mesh.nodes,N,beta);
 
 for i=1:N
        node_len=length(partmat(i).nodes);
        LHS_extra=sparse(node_len,node_len);
        for j=1:N
            if j==i
                continue
            else
                LHS_extra=LHS_extra+partmat(i).RHS_extra{j};
            end
        end
        partmat(i).LHS=alpha*LHS_extra;
        partmat(i).Gsub=partmat(i).Gsub+partmat(i).LHS;
 end

 tol = 1e-6;
maxit = 1000;

g_ini=cell(N,N); % Initial assumption for interface nodes

for i=1:N
    len_sub=length(partmat(i).nodes);
    for j=1:N
        if j==i
            continue
        else
            len_int=length(partmat(i).int_local{j});
            g_int=zeros(len_int,1);
            g_ini{i,j}=g_int;
            partmat(i).robin{j}=zeros(len_sub,1);
        end
    end
    
    partmat(i).uprev=zeros(len_sub,1);
    partmat(i).tol=1;
end

max_tol=1;
% max_tol=max([partmat(:).tol]);
tic
iter=0;
% while iter>0
 while max_tol>0.01 || iter<=2
    parfor i=1:N
        len_sub=length(partmat(i).nodes);
        RHS_extra=zeros(len_sub,1);
        for j=1:N
            if j==i
                continue
            else
               int_node=partmat(i).int_local{j};
               partmat(i).robin{j}(int_node)=g_ini{i,j};
            
               RHS_extra=RHS_extra+(partmat(i).RHS_extra{j}*partmat(i).robin{j});
            end
        end
        partmat(i).RHS_final=partmat(i).source+RHS_extra;
        partmat(i).soln=gmres(partmat(i).Gsub,partmat(i).RHS_final,[],tol,maxit);
%     partmat(i).soln=partmat(i).Gsub\partmat(i).RHS_final;
        partmat(i).tol=abs(partmat(i).soln-partmat(i).uprev)./abs(partmat(i).soln);
    end
    
   
    u_int=cell(N,N);
    for i=1:N
        for j=1:N
            if j==i
                continue
            else
                int_node=partmat(i).int_local{j};
                u_int{i,j}=partmat(i).soln(int_node);
            end
        end
    end
    
    g_next=cell(N,N);
    for i=1:N
        for j=1:N
            if j==i
                continue
            else 
                g_next{i,j}=2*alpha*u_int{j,i}-g_ini{j,i};
            end
        end
        partmat(i).uprev(:,:)=partmat(i).soln(:,:);
        partmat(i).maxtol=max(partmat(i).tol);
    end
    
    g_ini(:,:)=g_next(:,:);
  max_tol=max([partmat.maxtol]);
    iter=iter+1;
    disp(iter);
%      iter=iter-1;
 end
 
 toc
end