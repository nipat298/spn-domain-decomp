function Plot_mu2(P,T,h,muk,linestyle,cs)
%Function to plot output;
%P,T - FEM triangulation information
%muk- reconstructed parameter in elemental basis
%varargin = muaxf - actual parameter in elemental basis
%Limits for the domain
xmin=min(min(P(:,1)));
xmax=max(max(P(:,1)));
ymin=min(min(P(:,2)));
ymax=max(max(P(:,2)));
Nodes=size(P,1); Elem=size(T,1);

X=xmin:h:xmax;
Y=ymin:h:ymax;
Y=fliplr(Y);

xmu=zeros(Nodes,1);

for ii=1:Elem
    xmu(T(ii,:),1)=muk(1,1,ii);
end
for xx=1:length(X)
    for yy=1:length(Y)
        %             mu_rec2(xx,yy)=xmu((P(:,1)==X(xx))&(P(:,2)==Y(yy)));
        mu_rec(yy,xx)=xmu((P(:,1)==X(xx))&(P(:,2)==Y(yy)));
    end
end

fig2=figure();
imagesc(X,Y,(mu_rec))
%     set(gca,'CLim',[0,0.3],'FontName','TimesNewRoman','FontSize',13)
colorbar('FontName','TimesNewRoman','FontSize',13); colormap('jet')
axis('square')
xlabel('X co-ordinate in cms')
ylabel('Y co-ordinate in cms')



figure(200)
hold on
switch linestyle
    case 1
        plot(X,mu_rec(Y==cs,:),'-b')
    case 2
        plot(X,mu_rec(Y==cs,:),'-r')
    case 3
        plot(X,mu_rec(Y==cs,:),'-k')
    case 4
        plot(X,mu_rec(Y==cs,:),'--ko')
end
% figure(200)
% hold on
% plot(X,mu_rec(Y==-0.5,:))

% figure(200)
% hold on
% plot(X,mu_rec(29,:))

% figure(12)
% hold on
% plot(X,mu_rec(Y==0,:))


