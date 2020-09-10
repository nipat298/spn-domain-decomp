function S = GetSource(Dimension,N,Cgrad,CgradB,P,T,edgeElem,edgeType,Src)
%% Obtain the Source vector
%%%%%%%%%%%%%%% INPUT
% Dimension - Set 2 or 3 for 2D or 3D resp.
% N - order of SPn approximation
% Cgrad, CgradB - elemental coefficient matrices
% P - Node Coordinate matrix [NofNodes X 2] or [NofNodes X 3]
% T - Connectivity matrix. [NofElem X 3]
% edgeElem - [NofEdges X 1] - elements that lie on the domain edge
% edgeType - 1,2,3 depending on if the edge is side 1 side 2 or side 3 of
% the element
% Src - structure variable containing source information for NofSources

%%%%%%%%%%%%%%%%%%% OUTPUT
% S - [(N+1)/2 * NofNodes X NofSources] vector 
%% 
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%

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
    if strcmp(Src.Type,'beam')
        S=GetSourceX2D_beam(N,P,T,edgeElem,edgeType,Src,B);
    else
        if Dimension ==2
            S=GetSourceX2D(N,P,T,edgeElem,edgeType,Src,B);
        else
            S=GetSourceX3DN(N,P,T,edgeElem,edgeType,Src,B);
        end
    end
end