function [mesh,SrcData,Det,varargout]=MeshGen3D_new(saveFlag,Domain,h,SrcData,Det,varargin)
% Updated on 5th December 2020 - now supports extended sources.
% Modifications to Src and Det fields. Not backward compatible.
%Updated on 20 Nov 2020
%% Generating mesh data
xmin=Domain(1); xmax=Domain(2) ;ymin=Domain(3); ymax=Domain(4);zmin=Domain(5); zmax=Domain(6); % Defining the boundaries for the 2-D Phantom (cms)
mesh.hx=h(1); mesh.hy=h(2); mesh.hz = h(3);
[X,Y,Z]=meshgrid((xmin:mesh.hx:xmax),(ymin:mesh.hy:ymax),(zmin:mesh.hz:zmax));
mesh.nodes=[X(:),Y(:),Z(:)];          % x,y coordinates for the Nodes
mesh.tri=delaunayTriangulation(mesh.nodes); 
T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3),mesh.tri(:,4)];

%% Defining the inhomogeneity
if (nargin>5)
    mesh.inhom.elem=[];
    shape=varargin{1};
    centre=varargin{2};
    r=varargin{3};
    NofPhantom=size(centre,1);
    for ii=1:NofPhantom
        switch shape
            case 'cuboid'
                %% CUBOIDAL INHOMOGENEITY
                xc=centre(ii,1); yc=centre(ii,2); zc=centre(ii,3);
                len=r(ii,1); bre = r(ii,2); wid =r(ii,3);
                pF1= mesh.nodes((mesh.nodes(:,1)>=(xc-len/2))&(mesh.nodes(:,1)<=(xc+len/2))&(mesh.nodes(:,2)>=(yc-bre/2))&(mesh.nodes(:,2)<=(yc+bre/2))&(mesh.nodes(:,3)>=(zc-wid/2))&(mesh.nodes(:,3)<=(zc+wid/2)),:);
            case 'sphere'
                %% SPHERICAL INHOMOGENEITY
                xc=centre(ii,1); yc=centre(ii,2); zc=centre(ii,3);
                rad=r(ii,1); % radius of the circle
                dcircle=@(p) sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2 + (p(:,3)-zc).^2)-rad;
                pF1=mesh.nodes(feval(dcircle,mesh.nodes)<0,:); %Nodes on the fluorophore
        end
        for i=1:size(pF1,1)
            Nds_F1(i)=find((mesh.nodes(:,1)==pF1(i,1))&(mesh.nodes(:,2)==pF1(i,2))&(mesh.nodes(:,3)==pF1(i,3)));
        end
        % pF=P(Nds_F,:); %coordinates of points inside the fluorophore(inhomogeneity)
        BN=ismember(mesh.tri,Nds_F1);
        tF1=find(sum(BN,2)==4);
        mesh.inhom.ind(ii)=length(tF1);
        mesh.inhom.elem=[mesh.inhom.elem;tF1];
        clear Nds_F1 pF1 BN tF1
    end
end
%% Obtain Boundary edges
[mesh.edges,e_n]=boundedges_element_3D(mesh.nodes,T);
mesh.edgeElem=e_n(:,1);
mesh.edgeType=e_n(:,2);

%% Defining the source
NofSources =size(SrcData,2);

for s =1:NofSources
    if SrcData(s).ID=='S'
        % Surface/ Boundary source
        % find which surface/edge contains the source
        % presently supporting only circular/line sources sources
        % for simulating a point source, set source radius to less than
        % mesh.hx or mesh.hy
        xs = SrcData(s).Loc(1);ys = SrcData(s).Loc(2); zs = SrcData(s).Loc(3); % centre of the extended source/ location of point source
        if ~isfield(SrcData(s),'Radius')
            SrcData(s).Radius = 0.01*min(min(mesh.hx,mesh.hy),mesh.hz); % set default radius to smaller than mesh size . This simulates a point source
        end
        xyflag = (xs==xmin)||(xs==xmax); % '1' - source lies on x = constant edge ; '0' - source lies on y = constant edge or z = constant
        
        dcircle=@(p,r1,r2) sqrt((p(:,1)-r1).^2+(p(:,2)-r2).^2)-SrcData(s).Radius;
        if (xyflag)
            % Source lies on x= xmin // x = xmax
            SrcData(s).node = find((feval(dcircle,mesh.nodes(:,2:3),ys,zs)<0)&(mesh.nodes(:,1)==xs)); 
            % Node numbers
        else
             xyflag = (ys==ymin)||(ys==ymax);
            if (xyflag)
                % Source lies on y = ymin// y = ymax
                SrcData(s).node = find((feval(dcircle,mesh.nodes(:,[1,3]),xs,zs)<0)&(mesh.nodes(:,2)==ys)); 
            else
                % Source lies on z = ymin// z = ymax
                SrcData(s).node = find((feval(dcircle,mesh.nodes(:,1:2),xs,ys)<0)&(mesh.nodes(:,3)==zs)); 
            end
          
        end
        if isempty(SrcData(s).node)
            SrcData(s).node = nearestNeighbor(mesh.tri,xs,ys,zs);
        end
    else
        % Internal source
        % considers a circular source with a given radius.
        dcircle=@(p) sqrt((p(:,1)-xs).^2+(p(:,2)-ys).^2 + (p(:,3)-zs).^2)-SrcData(s).Radius;
        SrcData(s).node = find(feval(dcircle,mesh.nodes)<0);
    end
    SrcData(s).elem = mesh.edgeElem(sum(ismember(mesh.edges,SrcData(s).node),2)>0); % save all elements that contain the given node.
    % If a node belongs to two edges, both the edges / and hence both elems will be selected for that node.
    % This is different from previous versions, wherein either element was
    % selected.
end
%%
NofDet = size(Det.Loc,1);

if ~isfield(Det,'Radius')
    Det.Radius = 0.01*min(min(mesh.hx,mesh.hy),mesh.hz); % set default radius to simulate point source
end
xd = Det.Loc(:,1);yd = Det.Loc(:,2); zd = Det.Loc(:,3);
fprintf('Ensuring detectors are centred on Nodes... \n')
% push2node = nearestNeighbor(mesh.tri,xd,yd,zd);
% xd = mesh.nodes(push2node,1); yd = mesh.nodes(push2node,2);zd = mesh.nodes(push2node,3);

if Det.ID=='S'
    % Surface/ Boundary source
    % find which surface/edge contains the detector
    % presently supporting only circular/line sources sources
    % for simulating a point detector, set detector radius to less than
    % mesh.hx or mesh.hy
      dcircle=@(p,r1,r2) sqrt((p(:,1)-r1).^2+(p(:,2)-r2).^2)-Det.Radius;
    for d=1:NofDet
        xyflag = (xd(d)==xmin)||(xd(d)==xmax); % '1' - detector lies on x = constant edge ; '0' - detector lies on y = constant edge
      
        if (xyflag)
            % Source lies on x= xmin // x = xmax
            Det.node{d} = find((feval(dcircle,mesh.nodes(:,2:3),yd(d),zd(d))<0)&(mesh.nodes(:,1)==xd(d))); 
            
            % Node numbers
        else
             xyflag = (yd(d)==ymin)||(yd(d)==ymax);
            if (xyflag)
                % Source lies on y = ymin// y = ymax
                Det.node{d} = find((feval(dcircle,mesh.nodes(:,[1,3]),xd(d),zd(d))<0)&(mesh.nodes(:,2)==yd(d))); 
            else
                % Source lies on z = ymin// z = ymax
                Det.node{d} = find((feval(dcircle,mesh.nodes(:,1:2),xd(d),yd(d))<0)&(mesh.nodes(:,3)==zd(d))); 
            end
        end
        if isempty(Det.node{d})
            Det.node{d} =nearestNeighbor(mesh.tri,xd(d),yd(d),zd(d));
        end
        
                Det.elem{d} = mesh.edgeElem(sum(ismember(mesh.edges,Det.node{d}),2)>0);
    end
else
    % Internal source
    % considers a circular source with a given radius.
    for d=1:NofDet
        dcircle=@(p) sqrt((p(:,1)-xd(d)).^2+(p(:,2)-yd(d)).^2 + (p(:,3)-zd(d)).^2)-Det.Radius;
        Det.node{d} = find(feval(dcircle,mesh.nodes)<0);
        Det.elem{d} = mesh.edgeElem(sum(ismember(mesh.edges,Det.node{d}),2)>0);
    end
end


% remove detectors co-located with sources
% 
% ndet = length(Det.Loc);
% for s=1:NofSources
%     remove_det = zeros(ndet,1);
%     for d = 1:ndet
%     remove_det(d) = any(ismember(Det.elem{d},SrcData(s).elem));
%     end
%     ndet = ndet-sum(remove_det);
%     Det.elem = Det.elem(~remove_det);
%     Det.Loc = Det.Loc(~remove_det,:);
%     Det.node = Det.node(~remove_det);    
% end

%% Saving the mesh data
if (saveFlag)
clear xc yc r X Y dcircle Temp i ii r shape centre temp varargin
% fname=['mesh0_0',num2str(h*1e3),'.mat'];
fname=('mesh.mat');
save(fname,'mesh','SrcData','Det');
varargout{1}=fname;
end

