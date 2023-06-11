function [mesh,DT]=meshGen3D(Domain,mesh_spacing)

xmin=Domain(1); xmax=Domain(2);
ymin=Domain(3); ymax=Domain(4);
zmin=Domain(5); zmax=Domain(6);
[X,Y,Z] = meshgrid((xmin:mesh_spacing:xmax),(ymin:mesh_spacing:ymax),(zmin:mesh_spacing:zmax));
mesh.nodes= [X(:),Y(:),Z(:)];
DT =  delaunayTriangulation(mesh.nodes);
mesh.tri=[DT(:,1),DT(:,2),DT(:,3),DT(:,4)];
end
