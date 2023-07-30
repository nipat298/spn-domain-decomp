function t=boundtri_bdelem_3D(t,N)
%Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T) 
%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
%
%   Further modified and named as boundtri_bdelem_3D(). This function finds boundary triangles from tetrahedron mesh
%   and lists the element number corresponding to each boundary triangle

%   Inputs:
%   partmat: structure type containing subdomain info
%   NofDomains: Number of subdomains

%   Return variables:
%   partmat: Fields for boundary triangles and their corresponding elements added to 'partmat' structure  


% Form all edges, non-duplicates are boundary edges
%Form all possible triangles from the vertices/ nodes

for i=1:N
    edges=[t(i).ECM(:,[1,2,3]);
       t(i).ECM(:,[1,2,4]);
       t(i).ECM(:,[1,3,4]);
       t(i).ECM(:,[2,3,4])];
% List all vertices that have been excluded for the triangles above
node3=[t(i).ECM(:,4);t(i).ECM(:,3);t(i).ECM(:,2);t(i).ECM(:,1)];
% Sort all the edge triangles 
edges=sort(edges,2);
% Determine the unique triangles and their indices for the edges matrix
% stored in ix and for the foo matrix, stored in jx
[foo,ix,jx]=unique(edges,'rows','first');
% Find the number of occurrences of each triangle. Edge triangles will
% belong to a single tetrahedra, while unique triangles will lie on the
% edge
vec=histc(jx,1:max(jx));
qx=find(vec==1);
% retain only the triangles on the edge. 
e=edges(ix(qx),:); %also, e=foo(qx);
node3=node3(ix(qx)); %Determine the nodes that have been skipped
edge_no=ix(qx); %Store the edge number for the corresponding 
e_n=zeros(size(e,1),2);
i1=find(edge_no<=size(t(i).ECM,1));
e_n(i1,1)=edge_no(i1); 
e_n(i1,2)=1;
i2=find(edge_no(:)>size(t(i).ECM,1) & edge_no(:)<=2*size(t(i).ECM,1));
e_n(i2)=edge_no(i2)-size(t(i).ECM,1);
e_n(i2,2)=3;
i3=find(edge_no>2*size(t(i).ECM,1) & edge_no(:)<=3*size(t(i).ECM,1));
e_n(i3)=edge_no(i3)-2*size(t(i).ECM,1);
e_n(i3,2)=4;
i4=find(edge_no>3*size(t(i).ECM,1));
e_n(i4)=edge_no(i4)-3*size(t(i).ECM,1);
e_n(i4,2)=2;

t(i).boundtri=e;
t(i).belem=e_n;
end
end
