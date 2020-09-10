function Jmat = GetPartCurr(Dimension,N,CJ,P,T,edges,edgeElem,edgeType,Det)
if strcmp(Det.Type,'point')
    % point detector
    if Dimension ==2
        Jmat=GetExitCurrent2D(N,CJ,P,T,edgeElem,edgeType,Det);
    else
        Jmat=GetExitCurrent3D(N,CJ,P,T,edgeElem,edgeType,Det);
    end
else
    % extended detector  
    % Code not yet verified
%         if Dimension ==2
%             Jmat=GetExitCurr2D_extended(N,P,T,edges,edgeElem,edgeType,Det,CJ);
%         else
%             Jmat=GetExitCurr3D_extended(N,P,T,edges,edgeElem,edgeType,Det,CJ);
%         end
end