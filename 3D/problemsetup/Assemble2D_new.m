function fname=Assemble2D_new(mesh,varargin)
%Generate assembly matrices for Elemental Basis
%%%% INPUTS %%%%
% mesh - meshing file generated using MeshGen
%%%%%% OUTPUT %%%%%%%%%%
% fname - name of the mat file containing the assembly matrices. 
%%%%%%%%%%%%%%%%% REFERENCES %%%%%%%%%%%%%
% 1. The Finite Element Method in electromagnetics , Jiamming Jin
% 2. Omprakash Master's thesis
% 3. Coupled complex adjoint sensitivity, Fedele 2003
%%
% The Nodal basis is chosen as a+bx+cy;
P = mesh.nodes;
T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3)];
Element=size(T,1);

x1=P(T(:,1),1); y1=P(T(:,1),2);
x2=P(T(:,2),1); y2=P(T(:,2),2);
x3=P(T(:,3),1); y3=P(T(:,3),2);
%Generating the assembly matrices

b1(1,1,:)=y2-y3;           b2(1,1,:)=y3-y1;           b3(1,1,:)=y1-y2;
c1(1,1,:)=x3-x2;           c2(1,1,:)=x1-x3;           c3(1,1,:)=x2-x1;

Area=0.5*(b1.*c2-b2.*c1);
A=Area;

% stiffness matrix
Ks=[(b1.^2+c1.^2)./(4*A)     (b1.*b2+c1.*c2)./(4*A)   (b1.*b3+c1.*c3)./(4*A)
    (b2.*b1+c2.*c1)./(4*A)   (b2.^2+c2.^2)./(4*A)     (b2.*b3+c2.*c3)./(4*A)
    (b3.*b1+c3.*c1)./(4*A)   (b3.*b2+c3.*c2)./(4*A)   (b3.^2+c3.^2)./(4*A)];

% mass matrix
Km=[Area/6, Area/12, Area/12
    Area/12, Area/6, Area/12
    Area/12, Area/12, Area/6];

% [e,e_n]=boundedges_element_2D(P,T);
e=mesh.edges;
e_n(:,1)=mesh.edgeElem; e_n(:,2) = mesh.edgeType;
edge_l=zeros(1,1,Element); 
edge_l(1,1,e_n(:,1))=sqrt((P(e(:,1),1)-P(e(:,2),1)).^2+(P(e(:,1),2)-P(e(:,2),2)).^2);
% depending on edge type, assign appropriate boundary mass matrices to the
% corresponding edge element. In many works the boundary mass matrix is of
% size 2x2 for 2D, however we use 3x3 with the the quantities associated
% with the 3rd node set to zero for the corresponding edge depending on its
% edge type. This is for computational ease in the rest of the work, though
% this comes at the cost of additional memory.

f12=[1/3 1/6 0; 1/6 1/3 0; 0 0 0];
f13=[1/3 0 1/6; 0 0 0; 1/6 0 1/3];
f23=[0 0 0; 0 1/3 1/6; 0 1/6 1/3];


Kb=zeros(3,3,Element);
i1=find(e_n(:,2)==1);
Kb(:,:,e_n(i1,1))=repmat(f12,[1 1 size(i1,1)]).*repmat(edge_l(e_n(i1,1)),[3,3]);
i1=find(e_n(:,2)==2);
Kb(:,:,e_n(i1,1))=repmat(f23,[1 1 size(i1,1)]).*repmat(edge_l(e_n(i1,1)),[3,3]);
i1=find(e_n(:,2)==3);
Kb(:,:,e_n(i1,1))=repmat(f13,[1 1 size(i1,1)]).*repmat(edge_l(e_n(i1,1)),[3,3]);

% Iv,Jv are indexing matrices. These are helpful in global assembly as well
% as in the adjoint sensitivity computations. See Fedele2003 for details.

Iv=zeros(9,Element);Jv=zeros(9,Element);
for k=1:Element
    for j=1:3
    Iv(3*j-2:j*3,k)=T(k,:);
    Jv(3*j-2:j*3,k)=T(k,j);
    end
end
if isempty(varargin)
fname='assembly_data.mat';
else
    fname = varargin{1};
end
save(fname, 'Ks', 'Km', 'Kb', 'Iv','Jv','A','edge_l');

end
