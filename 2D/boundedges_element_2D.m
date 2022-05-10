function [e,e_n]=boundedges_element_2D(p,t)
%Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T) 
%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
%
%   Further modified and named as boundedges_element().
%   Return variables:
%   e: Boundary edges, represented by set of nodes number formed the edges.
%   e_n: Elements corresponding to boundary edges given by e array. 
%

% Form all edges, non-duplicates are boundary edges
edges=[t(:,[1,2]);
       t(:,[1,3]);
       t(:,[2,3])];
node3=[t(:,3);t(:,2);t(:,1)];
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows','first');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:); %also, e=foo(qx);
node3=node3(ix(qx));


edge_no=ix(qx);

e_n=zeros(size(e,1),2);
i1=find(edge_no<=size(t,1));
e_n(i1,1)=edge_no(i1); 
e_n(i1,2)=1;
i2=find(edge_no(:)>size(t,1) & edge_no(:)<=2*size(t,1));
e_n(i2)=edge_no(i2)-size(t,1);
e_n(i2,2)=3;
i3=find(edge_no>2*size(t,1));
e_n(i3)=edge_no(i3)-2*size(t,1);
e_n(i3,2)=2;


% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
iz=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(iz,[1,2])=e(iz,[2,1]);