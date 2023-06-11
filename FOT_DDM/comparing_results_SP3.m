% load('Phi values for single source');

%load('Phi_3D.mat');

load('Phi_nonDD_1.mat');
load('partmat_run1_tol0.005.mat');

u_em=PhiX;
chk1=abs(u_em(:));
figure()
plot(chk1(:));

% % EC=DT(:,:);
% % figure()
% % trisurf(EC,NodeCoord(:,1),NodeCoord(:,2),abs(u_DOT)),shading interp;
% % colorbar
% node_len=length(NodeCoord(:,1));
node_len=68921;
u_DD=zeros(2*node_len,1);
u_DD(partmat(1).nodes)=partmat(1).soln(1:length(partmat(1).nodes));
n11(:)=partmat(1).nodes(:)+ node_len;
n11=n11';
u_DD(n11)=partmat(1).soln((1+length(partmat(1).nodes)):end);

u_DD(partmat(2).nodes)=partmat(2).soln(1:length(partmat(2).nodes)); %Because nodes from second subdomain are being assigned value after first, so interface nodes would automatically take up the fluence values from 2nd subdomain
n22(:)=partmat(2).nodes(:)+ node_len;
n22=n22';
u_DD(n22)=partmat(2).soln((1+length(partmat(2).nodes)):end);

u_DD(partmat(3).nodes)=partmat(3).soln(1:length(partmat(3).nodes));
n33(:)=partmat(3).nodes(:)+ node_len;
n33=n33';
u_DD(n33)=partmat(3).soln((1+length(partmat(3).nodes)):end);

u_DD(partmat(4).nodes)=partmat(4).soln(1:length(partmat(4).nodes));
n44(:)=partmat(4).nodes(:)+ node_len;
n44=n44';
u_DD(n44)=partmat(4).soln((1+length(partmat(4).nodes)):end);
chk2=abs(u_DD(:));
figure()
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
