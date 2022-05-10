function fig2 = Plot_mu2(P,T,hx,hy,muk,linestyle,cs,varargin)
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

X=xmin:hx:xmax;
Y=ymin:hy:ymax;
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



figure(200)
hold on
switch linestyle
    case 1
        plot(X,mu_rec(Y==cs,:),'-b','LineWidth',1)
    case 2
        plot(X,mu_rec(Y==cs,:),'-r','LineWidth',1)
    case 3
        plot(X,mu_rec(Y==cs,:),'-k','LineWidth',1)
    case 4
        plot(X,mu_rec(Y==cs,:),'--ko','LineWidth',1)
end

if ~isempty(varargin)
    cs2 = varargin{1};
switch linestyle
    case 1
        p1 = plot(X,mu_rec(Y==cs2,:),'-b','LineWidth',1);
    case 2
        p1 = plot(X,mu_rec(Y==cs2,:),'-r','LineWidth',1);
    case 3
        p1 = plot(X,mu_rec(Y==cs2,:),'-k','LineWidth',1);
    case 4
        p1 = plot(X,mu_rec(Y==cs2,:),'--ko','LineWidth',1);
end
  p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
end


fig2=figure();
imagesc(X,Y,(mu_rec))
%     set(gca,'CLim',[0,0.3],'FontName','TimesNewRoman','FontSize',13)
colorbar('FontName','TimesNewRoman','FontSize',13); colormap('jet')
% axis('square')
xlabel('X co-ordinate in cms')
ylabel('Y co-ordinate in cms')



% figure(200)
% hold on
% plot(X,mu_rec(Y==-0.5,:))

% figure(200)
% hold on
% plot(X,mu_rec(29,:))

% figure(12)
% hold on
% plot(X,mu_rec(Y==0,:))


