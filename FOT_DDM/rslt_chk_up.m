% y1=load('Phi_nonDD_1');
 y1=load('PhiM_(-0.3,0,0).mat');
%y2=load('PhiX_DD_og2');
 u_em=y1.PhiM;
 u_DD=PhiM;
 
 rel_err=(abs(u_DD-u_em))./abs(u_em);
rel_err_cent=rel_err*100;
chk3=rel_err_cent(:);
figure()
plot(chk3(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 y1m=load('PhiM_run20_el3D');
 u_emm=y1m.PhiM;
 u_DDm=PhiM;
 
 
 y1=load('Phi_nonDD_1');
 y2=load('Phi_DD_1_tol0.002.mat');
 u_em=y1.PhiX;
 u_DD=y2.PhiX;
 
 Nodes=68921;
 
chk1=abs(u_em(1:Nodes,:));
figure()
plot(chk1(:));

hold on
chk2=abs(u_DD(1:Nodes,:));
% figure()
plot(chk2(:));

rel_err=(abs(u_DD-u_em))./abs(u_em);
rel_err_cent=rel_err*100;
chk3=rel_err_cent(:);
figure()
plot(chk3(:));

rel_err_m=(abs(u_DDm-u_emm))./abs(u_emm);
rel_err_cent_m=rel_err_m*100;
chk4=rel_err_cent_m(:);
figure()
plot(chk4(:));

tol1=1e-6;
 nn=find(abs(mesh.nodes(:,1)-(0.3))<tol1);
 
 nn=find(abs(mesh.nodes(:,3)==0));
 
 nn=find(abs(mesh.nodes(:,1)==0));
 
 nn=find(abs(mesh.nodes(:,2)==0));
 
 u_em_rr= reshape(u_em(nn),41,41);
 u_DD_rr=reshape(u_DD(nn),41,41);
 
 figure()
subplot(121)
imagesc(real(log(u_em_rr)));
colorbar
caxis([-7 7]);

subplot(122)
imagesc(real(log(u_DD_rr)));
colorbar
caxis([-7 7]);

chk4=abs(u_em(partmat(4).nodes));
figure()
plot(chk4(:));
hold on
chk5=abs(u_DD(partmat(4).nodes));
plot(chk5(:));



diff=u_em(:)-u_DD(:);
result=[u_em diff u_DD];


rel_err=(abs(u_DD(((1+68921):end),:)-u_em(((1+68921):end),:)))./abs(u_em(((1+68921):end),:));
rel_err_cent=rel_err*100;
chk3=rel_err_cent(:);
figure()
plot(chk3(:));

rel_err=(abs(u_DD((1:68921),:)-u_em((1:68921),:)))./abs(u_em((1:68921),:));
rel_err_cent=rel_err*100;
chk3=rel_err_cent(:);
figure()
plot(chk3(:));