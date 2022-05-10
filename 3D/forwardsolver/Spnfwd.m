function A=Spnfwd(Dimension,N,Cgrad,C,CgradB,CB,mesh,Iv,Jv,Ks,Km,Kb)
%%
% Subroutine to generate the FEM global assembly matrix from the
% coefficient matrices and the precomputed assemblymatrices.
%%
% load(AssemblyFile);
Element=size(mesh.tri,1);
Nodes=size(mesh.nodes,1);
B=zeros((N+1)/2,(N+1)/2,Element);    B1=B;
%Generates the boundary term in equation (57)
for i=1:size(mesh.edgeElem,1)
    B(:,:,mesh.edgeElem(i,1))=Cgrad(:,:,mesh.edgeElem(i,1))/CgradB(:,:,mesh.edgeElem(i,1));
    B1(:,:,mesh.edgeElem(i,1))=B(:,:,mesh.edgeElem(i,1))*CB(:,:,mesh.edgeElem(i,1));
end

 disp('Generate Sparse Assembly matrix... ')
     A = GetAssembledMat(Dimension,N,Iv,Jv,Nodes,Cgrad,Ks,C,Km,B1,Kb);
 
%  fprintf('Assembled matrix generated in %f seconds... \n',end_time-start_time);




