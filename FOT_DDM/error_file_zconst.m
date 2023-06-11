z1=1.8;
z2=1.9;

zz = [-1 1];
xy = [1 -1];

tol1=1e-6;
nn1=find((abs(mesh.nodes(:,3)-(z1))<tol1) & (abs(mesh.nodes(:,3)-(z1))>=0));
nn2=find((abs(mesh.nodes(:,3)-(z2))<tol1) & (abs(mesh.nodes(:,3)-(z2))>=0));

 u_rr1= reshape(rel_err_cent(nn1),41,41);
 u_rr2=reshape(rel_err_cent(nn2),41,41);
 
 figure()
subplot(121)
imagesc(zz,xy,real(log10(u_rr1)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along z=', num2str(z1), 'plane'])
xlabel('X')
ylabel('Y')
% axis([0 2 1 -1])
% c1 = caxis;


subplot(122)
imagesc(zz,xy,real(log10(u_rr2)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot along z=', num2str(z2), 'plane'])
xlabel('X')
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
title(['Percentage relative error plot(phi2) along z=', num2str(z1), 'plane'])
xlabel('X')
ylabel('Y')
% axis([0 2 1 -1])
% c1 = caxis;


subplot(122)
imagesc(zz,xy,real(log10(u_rr2_phi2)));
colorbar
caxis([-1 log10(max(rel_err_cent))]);
title(['Percentage relative error plot(phi2) along z=', num2str(z2), 'plane'])
xlabel('X')
ylabel('Y')
 