function A= AdjointFwd(Dimension,N,Cgrad,C,CgradB,CB,mesh,Ks,Km,Kb,Iv,Jv)
%%

Nodes=size(mesh.nodes,1); %No of nodes
Element=size(mesh.tri,1);
CBadj=zeros((N+1)/2,(N+1)/2,Element);
if N==1
    CBadj=permute(CB,[2,1,3]);
else
    for i=1:size(mesh.edgeElem,1)
        CBadj(:,:,mesh.edgeElem(i,1))=(Cgrad(:,:,mesh.edgeElem(i,1))*((CgradB(:,:,mesh.edgeElem(i,1)))\CB(:,:,mesh.edgeElem(i,1))));
    end
    CBadj=permute(CBadj,[2 1 3]);
    C=permute(C,[2 1 3]);
end

A=GetAssembledMat(Dimension,N,Iv,Jv,Nodes,conj(Cgrad),Ks,conj(C),Km,conj(CBadj),Kb);

end