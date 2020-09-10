function [mesh,SrcData,Det,varargout]=MeshGen2D_new(saveFlag,Domain,h,SrcData,Det,varargin)
%% Last updated on 29 November 2019
% Function to create triangulation mesh and obtain meshing and structural
% information
% Includes modification for considering source information as a structure
% variable
%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%
%saveFlag - '1' - save meshing data in mesh.mat file, else '0'
% Domain - [xmin,xmax,ymin,ymax] for the rectangular domain
% h - [hx,hy], mesh spacing
% see sample_el.m or sample_fl.m for details about the following input
% arguments
% SrcData - Structure containing source information.
% Det - Structure containing detector information.  
% varargin{1} - shape of the inclusion
% varagin{2} - centre of the inclusion
% varagin{3} - r, variable with supplementary information about the
% inclusion

%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh - structure containing meshing information
% SrcData - structure contianing source information
% Det - structure containing detector information
% if saveFlag = 1, varagout{1} contains name of the meshing file that is
% stored
%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
%% Generating mesh data
xmin=Domain(1); xmax=Domain(2) ;ymin=Domain(3); ymax=Domain(4); % Defining the boundaries for the 2-D Phantom (cms)
mesh.hx=h(1); mesh.hy=h(2); % mesh spacing
[X,Y]=meshgrid((xmin:mesh.hx:xmax),(ymin:mesh.hy:ymax));
mesh.nodes=[X(:),Y(:)];          % x,y coordinates for the Nodes
mesh.tri=delaunayTriangulation(mesh.nodes); % Triangulation class. see MATLAB's help in function for details
T=[mesh.tri(:,1),mesh.tri(:,2),mesh.tri(:,3)];

%% Defining the inhomogeneity
if nargin>5
    mesh.inhom.elem=[];
    shape=varargin{1};
    centre=varargin{2};
    r=varargin{3};
    NofPhantom=size(centre,1);
    for ii=1:NofPhantom
        switch shape
            case 'rect'
                %% SQUARE INHOMOGENEITY
                xc=centre(ii,1); yc=centre(ii,2);
                len=r(ii,1); bre = r(ii,2);
                pF1= mesh.nodes((mesh.nodes(:,1)>=(xc-len/2))&(mesh.nodes(:,1)<=(xc+len/2))&(mesh.nodes(:,2)>=(yc-bre/2))&(mesh.nodes(:,2)<=(yc+bre/2)),:);
            case 'ellipse'
                %% ELLIPTICAL INHOMOGENEITY
                xc=centre(ii,1);yc=centre(ii,2);
                major_axis=r(ii,1); minor_axis=r(ii,2);
                dellipse=@(p)((p(:,1)-xc)/major_axis).^2+((p(:,2)-yc)/minor_axis).^2-1;
                pF1=mesh.nodes(feval(dellipse,mesh.nodes)<0,:);
            case 'circle'
                %% CIRCULAR INHOMOGENEITY
                xc=centre(ii,1); yc=centre(ii,2);
                rad=r(ii,1); % radius of the circle
                dcircle=@(p) sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-rad;
                pF1=mesh.nodes(feval(dcircle,mesh.nodes)<0,:); %Nodes on the fluorophore
            case 'bean'
                xc=centre(ii,1); yc=centre(ii,2); %Center of the bean shaped inhomogeneity
                scale=r(ii,1);
                dbean=@(p) (p(:,1)-xc).^4+(p(:,2)-yc).^4 + ((p(:,1)-xc).^2).*((p(:,2)-yc).^2) - scale*(p(:,1)-xc).*((p(:,1)-xc).^2 + (p(:,2)-yc).^2);
                pF1=mesh.nodes(feval(dbean,mesh.nodes)<0,:); %Nodes on the fluorophore
                
        end
        for i=1:size(pF1,1)
            Nds_F1(i)=find((mesh.nodes(:,1)==pF1(i,1))&(mesh.nodes(:,2)==pF1(i,2)));
        end
        % pF=P(Nds_F,:); %coordinates of points inside the fluorophore(inhomogeneity)
        BN=ismember(mesh.tri,Nds_F1);
        tF1=find(sum(BN,2)==3);
        mesh.inhom.ind(ii)=length(tF1);
        mesh.inhom.elem=[mesh.inhom.elem;tF1];
        clear Nds_F1 pF1 BN tF1
    end
end
%% Defining the source
[mesh.edges,e_n]=boundedges_element_2D(mesh.nodes,T); %find edges and element numbers on the boundary
mesh.edgeElem=e_n(:,1);
mesh.edgeType=e_n(:,2);
NofSources =size(SrcData,2);
for s =1:NofSources
    if strcmp(SrcData(s).Type,'beam')
        % This section needs to be tested.
        eps =1e-12; % small push to ensure unique element is located.
        
        xc = SrcData(s).Loc(ii,1);yc = SrcData(s).Loc(ii,2);
        line_x = @(p) abs(p(:,1)-xc);
        line_y = @(p) abs(p(:,2)-yc);
        
        if SrcData(s).Radius(ii,1)==0
            pF=mesh.nodes((feval(line_y,mesh.nodes)<=SrcData(s).Radius(ii,2))&(mesh.nodes(:,1)==SrcData(s).Loc(ii,1)),:);
            querypoints =[pF(:,1),pF(:,2)+eps];
        else
            pF=mesh.nodes((feval(line_x,mesh.nodes)<=SrcData(s).Radius(ii,1))&(mesh.nodes(:,2)==SrcData(s).Loc(ii,2)),:);
            querypoints =[pF(:,1)+eps,pF(:,2)];
        end
        SrcData(s).beamelem(ii) = length(pF);
        
        SrcData(s).elem=pointLocation(mesh.tri,querypoints(:,1),querypoints(:,2));
        SrcData(s).node=pF;
        
    else
        if (SrcData(s).ID =='I')
            %internal source
            SrcData(s).elem = pointLocation(mesh.tri,SrcData(s).Loc);           
        else
            % Boundary source. check with loc in edge element 
            BaryCoords = zeros(length(mesh.edgeElem),3);
            x1 = mesh.nodes(mesh.tri(mesh.edgeElem,1),1); x2 = mesh.nodes(mesh.tri(mesh.edgeElem,2),1); x3 = mesh.nodes(mesh.tri(mesh.edgeElem,3),1);
            y1 = mesh.nodes(mesh.tri(mesh.edgeElem,1),2); y2 = mesh.nodes(mesh.tri(mesh.edgeElem,2),2); y3 = mesh.nodes(mesh.tri(mesh.edgeElem,3),2);
            det = (y2-y3).*(x1-x3) + (x3-x2).*(y1-y3);
            BaryCoords(:,1) = ((y2-y3).*(SrcData(s).Loc(1,1)-x3) + (x3-x2).*(SrcData(s).Loc(1,2)-y3))./det;
            BaryCoords(:,2) = ((y3-y1).*(SrcData(s).Loc(1,1)-x3) + (x1-x3).*(SrcData(s).Loc(1,2)-y3))./det;
            BaryCoords(:,3) = 1-sum(BaryCoords,2);
            BaryIn(:,1) = (BaryCoords(:,1)>=0)&(BaryCoords(:,1)<=1);
            BaryIn(:,2) = (BaryCoords(:,2)>=0)&(BaryCoords(:,2)<=1);
            BaryIn(:,3) = (BaryCoords(:,3)>=0)&(BaryCoords(:,3)<=1);
            ind = find(sum(BaryIn,2)>=3);
            SrcData(s).elem = mesh.edgeElem(ind(1));
        end
            SrcData(s).node = nearestNeighbor(mesh.tri,SrcData(s).Loc(:,1),SrcData(s).Loc(:,2));
    end
end
%% Defining the Detectors
% [e,e_n]=boundedges_element_2D(P,T);
NofDet = size(Det.Loc,1);
if strcmp(Det.Type,'point')
    for i=1:NofDet
        %     temp=find((abs(P(:,1)-DetLoc(i,1))==min(abs(P(:,1)-DetLoc(i,1))))&(abs(P(:,2)-DetLoc(i,2))==min(abs(P(:,2)-DetLoc(i,2)))));
        %     DNodes(i)=temp(1);
        
        if (Det.ID=='S')
            % detectors on surface
            Det.nodes(i)=nearestNeighbor(mesh.tri,Det.Loc(i,1),Det.Loc(i,2));
            dInd= mesh.edges(:,2)==Det.nodes(i);
            Det.elem(i)=mesh.edgeElem(dInd,1);
        else
            %internal detectors
            Det.nodes(i)=nearestNeighbor(mesh.tri,Det.Loc(i,1),Det.Loc(i,2));
            Det.elem(i) = pointLocation(mesh.tri,Det.Loc(i,:));
        end
        %
    end
    
else
    % For extended detectors,  the commented code below needs to be verified.
    %     dF1=[];
    %     Det.elem=[];
    %     for ii=1:NofDet
    %         xc = Det.Loc(ii,1);yc = Det.Loc(ii,2);
    %         line_x = @(p) abs(p(:,1)-xc);
    %         line_y = @(p) abs(p(:,2)-yc);
    %         % for boundary detectors
    %         %         boundflag = (xc==xmin)||(xc==xmax)||(yc==ymin)||(yc==ymax);
    %         if(Det.ID=='S')
    %             if Det.Radius(ii,1)==0
    %                 pF=find((feval(line_y,mesh.nodes)<=Det.Radius(ii,2))&(mesh.nodes(:,1)==Det.Loc(ii,1)));
    %                 dF1=[dF1;(pF)];
    %                 Det.beamelem(ii) = length(pF);
    %                 % find matching edge
    %                 detector_edge = ismember(pF,mesh.edges(:,2));
    %                 Det.elem=[Det.elem;mesh.edgeElem(detector_edge)];
    %
    %             else
    %                 pF=find((feval(line_x,mesh.nodes)<=Det.Radius(ii,1))&(mesh.nodes(:,2)==Det.Loc(ii,2)));
    %                 dF1=[dF1;(pF)];
    %                 Det.beamelem(ii) = length(pF);
    %                 % find matching edge
    %                 detector_edge = ismember(pF,mesh.edges(:,2));
    %                 Det.elem=[Det.elem;mesh.edgeElem(detector_edge)];
    %
    %             end
    %         else
    %             Det.elem=[];
    %             if Det.Radius(ii,1)==0
    %                 pF=find((feval(line_y,mesh.nodes)<=Det.Radius(ii,2))&(mesh.nodes(:,1)==Det.Loc(ii,1)));
    %                 dF1=[dF1;(pF)];
    %                 Det.beamelem(ii) = length(pF);
    %                 Det.elem=[Det.elem;pointLocation(mesh.tri,mesh.nodes(pF))];
    %             else
    %                 pF=find((feval(line_x,mesh.nodes)<=Det.Radius(ii,1))&(mesh.nodes(:,2)==Det.Loc(ii,2)));
    %                 dF1=[dF1;(pF)];
    %                 Det.beamelem(ii) = length(pF);
    %                 Det.elem=[Det.elem;pointLocation(mesh.tri,mesh.nodes(pF))];
    %             end
    %         end
    %     end
    %     Det.nodes=dF1;
end

% ignore_det_mask = zeros(size(Det.Loc,1),size(SrcData(s).elem,1));
% for i=1:NofSources
%     %     dist_frm_src = sqrt((Det.Loc(:,1)-SrcData.Loc(i,1)).^2+(Det.Loc(:,2)-SrcData.Loc(i,2)).^2);
%     dist_frm_src = sqrt((Det.Loc(:,1)-mesh.nodes(SrcData(s).node(i),1)).^2+(Det.Loc(:,2)-mesh.nodes(SrcData(s).node(i),2)).^2);
%     ignore_det_mask(:,i)= -2*(dist_frm_src<3*max(mesh.hx,mesh.hy))+1;
% end
% detE=pointLocation(DT,Det);
%% Saving the mesh data
if (saveFlag)
    clear xc yc r X Y dcircle Temp i ii r shape centre temp varargin
    % fname=['mesh0_0',num2str(h*1e3),'.mat'];
    fname=('mesh.mat');
    save(fname)
    varargout{1}=fname;
end

