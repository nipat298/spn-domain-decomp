function S = GetSource(Dimension,N,Cgrad,CgradB,P,T,Src)
% This function performs source calculation depending on whether the source is an internal source or surface source  

%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension: Choose 2 for 2D or 3 for 3D
% N: order of the SPN approximation (3 for SP3 case)
% Cgrad, CgradB: CgradphiX,CgradbphiX for the entire domain
% P: Nodal connectivity matrix for the entire domain
% T: Elemental Connectivity matrix for entire domain
% Src: source variable

%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S: Source vector

if Src.ID=='I'
    % Internal source
    
    if Dimension ==2
        S=GetSource2D_internal(N,P,T,Src);
    else
        S=GetSource3D_internal(N,P,T,Src);
    end
else
    % Boundary/ surface source
    B=zeros((N+1)/2,(N+1)/2,size(Src.elem,1));
    %Generates the boundary term in equation (57)
    for i=1:size(Src.elem,1)
        B(:,:,i)=Cgrad(:,:,Src.elem(i,1))/CgradB(:,:,Src.elem(i,1));
    end
        if Dimension ==2
%             S=GetSourceX2D(N,P,T,edgeElem,edgeType,Src,B);
             S=GetSourceX2D(N,P,T,Src,B);
        else
            S=GetSourceX3D(N,P,T,Src,B);
        end  
end