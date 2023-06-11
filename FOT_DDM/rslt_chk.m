 y1=load('PhiX_og2');
 y2=load('PhiX_DD_og2');
 u_em=y1.PhiX;
 u_DD=y2.PhiX;
 
chk1=abs(u_em(:));
figure()
plot(chk1(:));


chk2=abs(u_DD(:));
figure()
plot(chk2(:));


tol1=1e-6;
 nn=find(abs(mesh.nodes(:,1)-(-0.3))<tol1);
 
 nn=find(abs(mesh.nodes(:,3)==0));
 
 nn=find(abs(mesh.nodes(:,1)==0));
 
 u_em_rr= reshape(u_em(nn),41,41);
 u_DD_rr=reshape(u_DD(nn),41,41);
 
 figure()
subplot(121)
imagesc(real(log(u_em_rr)));

subplot(122)
imagesc(real(log(u_DD_rr)));

chk4=abs(u_em(partmat(4).nodes));
figure()
plot(chk4(:));
hold on
chk5=abs(u_DD(partmat(4).nodes));
plot(chk5(:));
