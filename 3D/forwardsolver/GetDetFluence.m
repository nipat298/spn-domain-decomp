function J = GetDetFluence(N,Nodes,Det,PhiX)

if strcmp(Det.Type,'point')
    switch (N)
        case 1
            J = PhiX(Det.Nodes);
        case 3
            J = PhiX(Det.Nodes)-(2/3)*PhiX(Det.Nodes+Nodes);
        case 5
            J = PhiX(Det.Nodes)-(2/3)*PhiX(Det.Nodes+Nodes) + (8/15)*PhiX(Det.Nodes+2*Nodes);
        case 7
            J = PhiX(Det.Nodes)-(2/3)*PhiX(Det.Nodes+Nodes) + (8/15)*PhiX(Det.Nodes+2*Nodes) - (16/35)*PhiX(Det.Nodes+3*Nodes);
    end
else
    
    elemcount=cumsum([0,Det.beamelem]);
    for ii=1:numel(elemcount)-1
        ind =elemcount(ii)+1:elemcount(ii);
        switch (N)
            case 1
                J(ii) = sum(PhiX(Det.Nodes(ind)));
            case 3
                J(ii) = sum(PhiX(Det.Nodes(ind))-(2/3)*PhiX(Det.Nodes(ind)+Nodes));
            case 5
                J = sum(PhiX(Det.Nodes(ind))-(2/3)*PhiX(Det.Nodes(ind)+Nodes) + (8/15)*PhiX(Det.Nodes(ind)+2*Nodes));
            case 7
                J = sum(PhiX(Det.Nodes(ind))-(2/3)*PhiX(Det.Nodes(ind)+Nodes) + (8/15)*PhiX(Det.Nodes(ind)+2*Nodes) - (16/35)*PhiX(Det.Nodes(ind)+3*Nodes));
        end
    end
end