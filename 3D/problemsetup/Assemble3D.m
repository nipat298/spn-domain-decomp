function fname=Assemble3D(mesh,varargin)
%Generate assembly matrices for Elemental Basis
% The Nodal basis is chosen as a+bx+cy+dz;
% 

P = mesh.nodes;
T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3),mesh.tri(:,4)];
Element=size(T,1);

x1=P(T(:,1),1); y1=P(T(:,1),2); z1=P(T(:,1),3);
x2=P(T(:,2),1); y2=P(T(:,2),2); z2=P(T(:,2),3);
x3=P(T(:,3),1); y3=P(T(:,3),2); z3=P(T(:,3),3);
x4=P(T(:,4),1); y4=P(T(:,4),2); z4=P(T(:,4),3);
%Generating the assembly matrices
% V=zeros(1,1,Element);
% for ii=1:Element
% V(1,1,ii)=(1/6)*det([1,1,1,1; x1(ii),x2(ii),x3(ii),x4(ii); y1(ii),y2(ii),y3(ii),y4(ii);z1(ii),z2(ii),z3(ii),z4(ii)]);
% end
% V=abs(V);

s1 = y3.*z4 - y4.*z3; s2 = y2.*z4 - y4.*z2; s3 = y2.*z3 - y3.*z2; 
s4 = y1.*z4 - y4.*z1; s6 = y1.*z2 - y2.*z1; s5 = y1.*z3 - y3.*z1;

V(1,1,:) = (1/6)*(x2.*s1 - x3.*s2 + x4.*s3...
    - x1.*s1 + x3.*s4 - x4.*s5...
    + x1.*s2 - x2.*s4 + x4.*s6...
    - x1.*s3 + x2.*s5 - x3.*s6);



b1(1,1,:)=-(s1 - s2 + s3);
b2(1,1,:)=s1 - s4 + s5;
b3(1,1,:)=-(s2 - s4 + s6);
b4(1,1,:)=s3 - s5 + s6;  

clear s1 s2 s3 s4 s5 s6

c1(1,1,:)=(x3.*z4-x4.*z3) - (x2.*z4-x4.*z2) + (x2.*z3-x3.*z2);
c2(1,1,:)=-((x3.*z4-x4.*z3) - (x1.*z4-x4.*z1) + (x1.*z3-x3.*z1));
c3(1,1,:)=(x2.*z4-x4.*z2) - (x1.*z4-x4.*z1) + (x1.*z2-x2.*z1);
c4(1,1,:)=-((x2.*z3-x3.*z2) - (x1.*z3-x3.*z1) + (x1.*z2-x2.*z1));  

d1(1,1,:)=-((x3.*y4-x4.*y3) - (x2.*y4-x4.*y2) + (x2.*y3-x3.*y2));
d2(1,1,:)=((x3.*y4-x4.*y3) - (x1.*y4-x4.*y1) + (x1.*y3-x3.*y1));
d3(1,1,:)=-((x2.*y4-x4.*y2) - (x1.*y4-x4.*y1) + (x1.*y2-x2.*y1));
d4(1,1,:)=((x2.*y3-x3.*y2) - (x1.*y3-x3.*y1) + (x1.*y2-x2.*y1));  

% Ks=[(b1.^2+c1.^2 + d1.^2)./(6*V)     (b1.*b2+c1.*c2 + d1.*d2)./(6*V)    (b1.*b3+c1.*c3 + d1.*d3)./(6*V) (b1.*b4 + c1.*c4 + d1.*d4)./(6*V)
%     (b2.*b1+c2.*c1 + d1.*d2)./(6*V)  (b2.^2+c2.^2 + d2.^2)./(6*V)       (b2.*b3+c2.*c3 + d2.*d3)./(6*V) (b2.*b4 + c2.*c4 + d2.*d4)./(6*V)
%     (b3.*b1+c3.*c1 + d1.*d3)./(6*V)  (b3.*b2+c3.*c2+d2.*d3)./(6*V)      (b3.^2+c3.^2 + d3.^2)./(6*V)    (b3.*b4 + c3.*c4 + d3.*d4)./(6*V)
%     (b4.*b1+c4.*c1 + d1.*d4)./(6*V)  (b4.*b2+c4.*c2+d2.*d4)./(6*V)      (b3.*b4+c3.*c4+d3.*d4)./(6*V)   (b4.^2+c4.^2 + d4.^2)./(6*V)];
% 
Ks=[(b1.^2+c1.^2 + d1.^2)./(36*V)     (b1.*b2+c1.*c2 + d1.*d2)./(36*V)    (b1.*b3+c1.*c3 + d1.*d3)./(36*V) (b1.*b4 + c1.*c4 + d1.*d4)./(36*V)
    (b2.*b1+c2.*c1 + d1.*d2)./(36*V)  (b2.^2+c2.^2 + d2.^2)./(36*V)       (b2.*b3+c2.*c3 + d2.*d3)./(36*V) (b2.*b4 + c2.*c4 + d2.*d4)./(36*V)
    (b3.*b1+c3.*c1 + d1.*d3)./(36*V)  (b3.*b2+c3.*c2+d2.*d3)./(36*V)      (b3.^2+c3.^2 + d3.^2)./(36*V)    (b3.*b4 + c3.*c4 + d3.*d4)./(36*V)
    (b4.*b1+c4.*c1 + d1.*d4)./(36*V)  (b4.*b2+c4.*c2+d2.*d4)./(36*V)      (b3.*b4+c3.*c4+d3.*d4)./(36*V)   (b4.^2+c4.^2 + d4.^2)./(36*V)];



Km=[V/10 V/20 V/20 V/20
    V/20 V/10 V/20 V/20
    V/20 V/20 V/10 V/20
    V/20 V/20 V/20 V/10];

e = mesh.edges; 
e_n(:,1) = mesh.edgeElem; e_n(:,2) = mesh.edgeType;
% [e,e_n]=boundedges_element_3D(P,T);
edge_A=zeros(1,1,Element); 
x1b=(P(e(:,1),1)); y1b=(P(e(:,1),2)); z1b=(P(e(:,1),3));
x2b=(P(e(:,2),1)); y2b=(P(e(:,2),2)); z2b=(P(e(:,2),3));
x3b=(P(e(:,3),1)); y3b=(P(e(:,3),2)); z3b=(P(e(:,3),3));

edge_A(1,1,e_n(:,1))=0.5.*sqrt(((y2b-y1b).*(z3b-z1b)-(z2b-z1b).*(y3b-y1b)).^2 + ((x2b-x1b).*(z3b-z1b)-(z2b-z1b).*(x3b-x1b)).^2 + ((x2b-x1b).*(y3b-y1b)-(y2b-y1b).*(x3b-x1b)).^2);


f123=[1/6 1/12 1/12 0; 1/12 1/6 1/12 0;1/12 1/12 1/6 0;0 0 0 0];
f124=[1/6 1/12 0 1/12; 1/12 1/6 0 1/12; 0 0 0 0; 1/12 1/12 0 1/6];
f134=[1/6 0 1/12 1/12; 0 0 0 0; 1/12 0 1/6 1/12; 1/12 0 1/12 1/6];
f234=[0 0 0 0; 0 1/6 1/12 1/12; 0 1/12 1/6 1/12; 0 1/12 1/12 1/6];

Kb=zeros(4,4,Element);
i1=find(e_n(:,2)==1);
Kb(:,:,e_n(i1,1))=repmat(f123,[1 1 size(i1,1)]).*repmat(edge_A(e_n(i1,1)),[4,4]);
i1=find(e_n(:,2)==2);
Kb(:,:,e_n(i1,1))=repmat(f234,[1 1 size(i1,1)]).*repmat(edge_A(e_n(i1,1)),[4,4]);
i1=find(e_n(:,2)==3);
Kb(:,:,e_n(i1,1))=repmat(f124,[1 1 size(i1,1)]).*repmat(edge_A(e_n(i1,1)),[4,4]);
i1=find(e_n(:,2)==4);
Kb(:,:,e_n(i1,1))=repmat(f134,[1 1 size(i1,1)]).*repmat(edge_A(e_n(i1,1)),[4,4]);

Iv=zeros(16,Element);Jv=zeros(16,Element);
for k=1:Element
    for j=1:4
    Iv(4*j-3:j*4,k)=T(k,:); Jv(4*j-3:j*4,k)=T(k,j);
    end
end

fname='assembly_data.mat';
save(fname, 'Ks', 'Km', 'Kb', 'Iv','Jv','V','edge_A','-v7.3');

end
