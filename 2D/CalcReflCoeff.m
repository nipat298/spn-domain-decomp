function [A,B,C,D,varargout]=CalcReflCoeff(N,n0,ni)
%% Function to calculate the reflection coefficient dependent quantities as
% defined in Klose 2005. 
%%%%%%%%%%%% INPUT
% N - order of SPN approximation
% n0 - refractive index of outside medium / fiber
% ni - refractive index of the medium
%%%%%%%%%%%% OUTPUT
% coefficients A,B,C,D, E, F, G, H, of size [(N+1)/2 X 1] each. 
% Refer Klose et al for details. 

%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com

%%
for ii=1:2*N
a=getr(ii,n0,ni);                                                
R(ii)=(integral(a,0,pi/2));
end

%Reflection parameters. ref. Klose 
A(1)=-R(1);
B(1)=3*R(2);

if N>1
    
    A(2)=(-2.25*R(1))+(7.5*R(3))-(6.25*R(5));
B(2)=(15.75*R(2))-(52.5*R(4))+(43.75*R(6));

C(1)=(-1.5*R(1))+(2.5*R(3));
D(1)=(1.5*R(2))-(2.5*R(4));


C(2)=C(1);
D(2)=D(1);
end

if N>3
    A(3)=-(225/64)*R(1) +(525/16)*R(3) -(3395/32)*R(5)+(2205/16)*R(7)-(3969/64)*R(9);
    B(3)=(2475/64)*R(2)-(5775/16)*R(4)+(37345/32)*R(6)-(24255/16)*R(8)+(43659/64)*R(10);
    C(3)=(15/8)*R(1)-(35/4)*R(3)+(63/8)*R(5);
    D(3)=-(15/8)*R(2)+(35/4)*R(4)-(63/8)*R(6);
    E(3)=-(45/16)*R(1) + (285/16)*R(3) - (539/16)*R(5) + (315/16)*R(7);
    F(3)=(45/16)*R(2)- (285/16)*R(4) + (539/16)*R(6)- (315/16)*R(8);
  
    E(1)=C(3);
    F(1)=D(3);
    E(2)=E(3);
    F(2)=F(3);   
    
    varargout{1}=E;
    varargout{2}=F;
end
if N>5
    G(3)=-(525/128)*R(1)+(7175/128)*R(3)-(17325/64)*R(5)+(37395/64)*R(7)-(73689/128)*R(9)+(27027/128)*R(11);
    H(3)=(525/128)*R(2)-(7175/128)*R(4)+(17325/64)*R(6)-(37395/64)*R(8)+(73689/128)*R(10)-(27027/128)*R(12);
    
    A(4)=-(1225/256)*R(1)+(11025/128)*R(3)-(147735/256)*R(5)+(116655/64)*R(7)-(750519/256)*R(9)+(297297/128)*R(11)-(184041/256)*R(13);
    B(4)=15*((1225/256)*R(2)-(11025/128)*R(4)+(147735/256)*R(6)-(116655/64)*R(8)+(750519/256)*R(10)-(297297/128)*R(12)+(184041/256)*R(14));
    C(4)=-(35/16)*R(1)+(315/16)*R(3)-(693/16)*R(5)+(429/16)*R(7);
    D(4)=(35/16)*R(2)-(315/16)*R(4)+(693/16)*R(6)-(429/16)*R(8);  
    E(4)=105/32*R(1)-35*R(3) + (1827/16)*R(5)-(297/2)*R(7)+(2145/32)*R(9);
    F(4)=-105/32*R(2)+35*R(4) - (1827/16)*R(6)+(297/2)*R(8)-(2145/32)*R(10);
    G(4)=G(3);
    H(4)=H(3);
    G(1)=C(4);
    H(1)=D(4);
    G(2)=E(4);
    H(2)=F(4);   
    
    varargout{1}=E;
    varargout{2}=F;
    varargout{3}=G;
    varargout{4}=H;

end

end

function rcoeff=getr(n,n0,ni)
%function returns a function handle to evaluate the 'n'th moment of the reflection coefficient for
%given n0 and ni
rcoeff=@refcffni;

    function Rn=refcffni(v1)
    vc=asin(n0/ni);
    temp=ones(size(v1));
    %Eqn (5) Klose_2006
    %for defns. of variables chk Klose paper
    ind=find(abs(v1)<vc);
    sinv2=(ni/n0)*sin(v1(ind));
    cosv2=sqrt(1-(sinv2).^2);
    temp(ind)=0.5*(((ni*cosv2-n0*cos(v1(ind)))./(ni*cosv2+n0*cos(v1(ind)))).^2)+...
        0.5*(((ni*cos(v1(ind))-n0*cosv2)./(ni*cos(v1(ind))+n0*cosv2)).^2);
   Rn=(temp).*(((cos(v1)).^n).*sin(v1));

    end
end