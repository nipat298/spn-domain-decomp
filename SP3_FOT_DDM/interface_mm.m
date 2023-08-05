function ABC=interface_mm(partmat,NodeCoord,N,beta)

% This function gives the assembled matrix corresponding to interface nodes for each subdomain

%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partmat: structure type variable containing subdomain info
% NodeCoord: nodal coordinate matrix for the entire undivided domain 
% N: Number of subdomains
% beta: penalty coeffcient (extra parameter along the interface) defined by Deng in his approach (method 2)

%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
% ABC: Structure type variable with fields corresponding to interface matrix for each interface of a given subdomain
%       (Example: Subdomain 2 would have interfaces along j=1 & 3. Therefore, 1st and 3rd row would be filled for subdomain 2.
%          Subdomain 1 would have interface only along j=1. Therefore, only 1st row would be filled for subdomain 1) 

%%% Variables in this function:
% elem_len: number of elements in that specific subdomain
% tri_int_area: stores the value of area of the triangles lying on the considered interface
% partmat.artbound: lists all the interface triangles(and their corresponding nodes) for the given subdomain
% partmat.art_tri_SD{j}: lists out the triangles(and their corresponding nodes)lying on the artificial boundary "ij"
% partmat.art_elem_SD{j}: lists out the elements correspoding to the triangles lying on the artificial boundary "ij"
% bound_tri(i): structure formed to store the 4 values required for evaluating area of a triangle
%%% The four values are a,b,c and s. 


NC=length(NodeCoord(:,1));
 elem_len=length(partmat.ECM(:,1));
     for j=1:N
            tri_int_area=zeros(1,1,elem_len); 
            ismem=ismember(partmat.art_bound,partmat.interface_g{j});
            indx=(ismem(:,1)==1 & ismem(:,2)==1 & ismem(:,3)==1);
            partmat.art_tri_SD{j}=partmat.art_bound(indx,:);
            partmat.art_elem_SD{j}=partmat.artb_elem(indx,:);
            bound_tri.int_a{j}=sqrt((NodeCoord(partmat.art_tri_SD{j}(:,1),1)-NodeCoord(partmat.art_tri_SD{j}(:,2),1)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,1),2)-NodeCoord(partmat.art_tri_SD{j}(:,2),2)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,1),3)-NodeCoord(partmat.art_tri_SD{j}(:,2),3)).^2);
            bound_tri.int_b{j}=sqrt((NodeCoord(partmat.art_tri_SD{j}(:,2),1)-NodeCoord(partmat.art_tri_SD{j}(:,3),1)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,2),2)-NodeCoord(partmat.art_tri_SD{j}(:,3),2)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,2),3)-NodeCoord(partmat.art_tri_SD{j}(:,3),3)).^2);
            bound_tri.int_c{j}=sqrt((NodeCoord(partmat.art_tri_SD{j}(:,1),1)-NodeCoord(partmat.art_tri_SD{j}(:,3),1)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,1),2)-NodeCoord(partmat.art_tri_SD{j}(:,3),2)).^2+(NodeCoord(partmat.art_tri_SD{j}(:,1),3)-NodeCoord(partmat.art_tri_SD{j}(:,3),3)).^2);
            bound_tri.int_s{j}=(bound_tri.int_a{j}+bound_tri.int_b{j}+bound_tri.int_c{j})/2;

            
            %%% Here we apply the formula for the area of a trianlge in terms of the length of its sides
            %%% area of triangle with sides a,b,c is : sqrt(s*(s-a)*(s-b)*(s-c))  
            tri_int_area(1,1,partmat.art_elem_SD{j}(:,1))=sqrt((bound_tri.int_s{j}).*(bound_tri.int_s{j}-bound_tri.int_a{j}).*(bound_tri.int_s{j}-bound_tri.int_b{j}).*(bound_tri.int_s{j}-bound_tri.int_c{j}));

            f_123=[1/6 1/12 1/12 0; 1/12 1/6 1/12 0; 1/12 1/12 1/6 0; 0 0 0 0];
            f_124=[1/6 1/12 0 1/12; 1/12 1/6 0 1/12; 0 0 0 0; 1/12 1/12 0 1/6];
            f_134=[1/6 0 1/12 1/12; 0 0 0 0; 1/12 0 1/6 1/12; 1/12 0 1/12 1/6];
            f_234=[0 0 0 0; 0 1/6 1/12 1/12; 0 1/12 1/6 1/12; 0 1/12 1/12 1/6];

            Kb_int=zeros(4,4,elem_len);

            ind=find(partmat.art_elem_SD{j}(:,2)==1);
            Kb_int(:,:,partmat.art_elem_SD{j}(ind,1))=repmat(f_123,[1 1 size(ind,1)]).*repmat(tri_int_area(partmat.art_elem_SD{j}(ind,1)),[4,4]);
            ind=find(partmat.art_elem_SD{j}(:,2)==2);
            Kb_int(:,:,partmat.art_elem_SD{j}(ind,1))=repmat(f_234,[1 1 size(ind,1)]).*repmat(tri_int_area(partmat.art_elem_SD{j}(ind,1)),[4,4]);
            ind=find(partmat.art_elem_SD{j}(:,2)==3);
            Kb_int(:,:,partmat.art_elem_SD{j}(ind,1))=repmat(f_124,[1 1 size(ind,1)]).*repmat(tri_int_area(partmat.art_elem_SD{j}(ind,1)),[4,4]);
            ind=find(partmat.art_elem_SD{j}(:,2)==4);
            Kb_int(:,:,partmat.art_elem_SD{j}(ind,1))=repmat(f_134,[1 1 size(ind,1)]).*repmat(tri_int_area(partmat.art_elem_SD{j}(ind,1)),[4,4]);

                       
            Ti = reshape((partmat.ECM(:,1:4))',[4,1,elem_len]);
            Tj = reshape((partmat.ECM(:,1:4))',[1,4,elem_len]);
            I = repmat(Ti,1,4);
            J = repmat(Tj,4,1);

            added_term = sparse(I(:),J(:),Kb_int(:),NC,NC);
            
            nodes_tot=1:length(NodeCoord);
            sub_ind=setdiff((nodes_tot)',partmat.nodes);
            
            added_term(sub_ind,:)=[];
            added_term(:,sub_ind)=[];
            
            RHS_extra=(1/beta)*added_term;
            ABC.RHS_extra{j}=RHS_extra;
    end
    
end


