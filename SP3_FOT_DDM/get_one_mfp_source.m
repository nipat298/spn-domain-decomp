function Src = get_one_mfp_source(Dimension,mesh,Src,mfp)

%Updated on 4 Decemeber 2020
% Now supports point and extended sources

if (Dimension==2)
fprintf('Moving source to 1mfp inside the medium.... \n');
ang = Src.ang_in*pi/180;
 if (mod(Src.ang_in*pi/180,pi/2)==0)&&(mod(Src.ang_in*pi/180,pi)~=0)
        dxs= cos(ang);
        dys = sin(ang);
    else
        dxs = sin(ang);
        dys = cos(ang);
    end

Src.Loc = mesh.nodes(Src.node,:) + [mfp*dxs,mfp*dys]; % update source location to one mfp inside the medium
Src.node = nearestNeighbor(mesh.tri,Src.Loc(:,1),Src.Loc(:,2));
Src.elem = pointLocation(mesh.tri,Src.Loc(:,1),Src.Loc(:,2));
Src.Strength = mean(mesh.opt.musx*(1-mesh.opt.gx))*mfp; % Added on 31102020 - this scaling is from Farrell+ Patterson (Med Phys. 1992)

else

%Updated on 4 Decemeber 2020
% Now supports point and extended sources
fprintf('Moving source to 1mfp inside the medium.... \n');

theta = Src.ang_in*pi/180;
phi = 0; % Presently phi is set to 0. should be able to extend support for arbitrary phi in future

% Only supports sources in y = constant. Extend to other planes in future
 if ((mod(theta,pi/2)==0)&&(theta~=0))
        dxs= cos(theta);
        dzs = sin(theta)*cos(phi);
        dys = sin(theta)*sin(phi);
        
  else
        dxs = sin(theta)*cos(phi);
        dys = sin(theta)*sin(phi);
        dzs = cos(theta);
  end



Src.Loc = mesh.nodes(Src.node,:) + mfp*[dxs,dys,dzs]; % update source location to one mfp inside the medium
Src.node = nearestNeighbor(mesh.tri,Src.Loc(:,1),Src.Loc(:,2),Src.Loc(:,3));
Src.elem = pointLocation(mesh.tri,Src.Loc(:,1),Src.Loc(:,2),Src.Loc(:,3));
Src.Strength = mean(mesh.opt.musx*(1-mesh.opt.gx))*mfp; % Added on 31102020 - this scaling is from Farrell+ Patterson (Med Phys. 1992)

end