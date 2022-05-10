function A=Spnfwd(Dimension,N,Cgrad,C,CgradB,CB,P,T,edgeElem,AssemblyFile)
%%
% Subroutine to generate the FEM global assembly matrix from the
% coefficient matrices and the precomputed assemblymatrices.

%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
load(AssemblyFile);
Element=size(T,1);
Nodes=size(P,1);
B=zeros((N+1)/2,(N+1)/2,Element);    B1=B;
%Generates the boundary term in equation (57)
for i=1:size(edgeElem,1)
    B(:,:,edgeElem(i,1))=Cgrad(:,:,edgeElem(i,1))/CgradB(:,:,edgeElem(i,1));
    B1(:,:,edgeElem(i,1))=B(:,:,edgeElem(i,1))*CB(:,:,edgeElem(i,1));
end
A=GetAssembledMat(Dimension,N,Iv,Jv,Nodes,Cgrad,Ks,C,Km,B1,Kb);

end


