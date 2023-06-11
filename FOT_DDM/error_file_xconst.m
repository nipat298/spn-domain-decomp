x1=0.2;
x2=0.8;

zz = [0 2];
xy = [1 -1];

tol1=1e-6;
nn1=find((abs(mesh.nodes(:,1)-(x1))<tol1) & (abs(mesh.nodes(:,1)-(x1))>=0));
nn2=find((abs(mesh.nodes(:,1)-(x2))<tol1) & (abs(mesh.nodes(:,1)-(x2))>=0));

 u_rr1= reshape(rel_err_cent(nn1),41,41);
 u_rr2=reshape(rel_err_cent(nn2),41,41);
 
 figure()
subplot(121)
imagesc(zz,xy,real(log10(u_rr1)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along x=',num2str(x1), 'plane'])
xlabel('Z')
ylabel('Y')
% axis([0 2 1 -1])
% c1 = caxis;


subplot(122)
imagesc(zz,xy,real(log10(u_rr2)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along x=',num2str(x2), 'plane'])
xlabel('Z')
ylabel('Y')


%%%%Plotting for phi2

nn1_phi2=nn1+Nodes;
nn2_phi2=nn2+Nodes;

 u_rr1_phi2= reshape(rel_err_cent(nn1_phi2),41,41);
 u_rr2_phi2=reshape(rel_err_cent(nn2_phi2),41,41);

figure()
subplot(121)
imagesc(zz,xy,real(log10(u_rr1_phi2)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along x=',num2str(x1), 'plane'])
xlabel('Z')
ylabel('Y')
% axis([0 2 1 -1])
% c1 = caxis;


subplot(122)
imagesc(zz,xy,real(log10(u_rr2_phi2)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along x=',num2str(x2), 'plane'])
xlabel('Z')
ylabel('Y')

