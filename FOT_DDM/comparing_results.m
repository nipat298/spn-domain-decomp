% load('Phi values for single source');

%load('Phi_3D.mat');

y1=load('PhiX.mat');
y2=load('PhiX_DD.mat');

u_em=y1.PhiX;
u_DD=y2.PhiX;

chk1=abs(u_em(:));
figure()
plot(chk1(:));

% % EC=DT(:,:);
% % figure()
% % trisurf(EC,NodeCoord(:,1),NodeCoord(:,2),abs(u_DOT)),shading interp;
% % colorbar
% node_len=length(NodeCoord(:,1));
% node_len=68921;
% u_DD=zeros(node_len,1);
% 
% u_DD(partmat(1).nodes)=partmat(1).soln;
% u_DD(partmat(2).nodes)=partmat(2).soln; %Because nodes from second subdomain are being assigned value after first, so interface nodes would automatically take up the fluence values from 2nd subdomain
% u_DD(partmat(3).nodes)=partmat(3).soln;
% u_DD(partmat(4).nodes)=partmat(4).soln;

hold on
chk2=abs(u_DD(:));

plot(chk2(:));

% u_DD(value00)=u0;
% u_DD(value11)=u1; %Because nodes from second subdomain are being assigned value after first, so interface nodes would automatically take up the fluence values from 2nd subdomain
% 

diff=u_em(:)-u_DD(:);

result=[u_em diff u_DD];
% % 
% % figure()
% % trisurf(EC,NodeCoord(:,1),NodeCoord(:,2),abs(u_DD)),shading interp;
% % colorbar
hold on
rel_err=(abs(u_DD-u_em))./abs(u_em);
rel_err_cent=rel_err*100;
chk3=rel_err_cent(:);
figure()
plot(chk3(:));

diff_per=abs(diff);
chk4=diff(:);
figure()
plot(chk4(:));
% % figure()
% % trisurf(EC,NodeCoord(:,1),NodeCoord(:,2),rel_err*100),shading interp;
% % colorbar
% % 
% % chk=real(u_DOT);
% % chk2=real(u_DD);

% % figure
% % plot(chk./chk2);

chk4=abs(u_em(partmat(4).nodes));
figure()
plot(chk4(:));

hold on



chk5=abs(u_DD(partmat(4).nodes));
plot(chk5(:));

