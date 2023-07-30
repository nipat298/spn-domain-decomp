% function [A,B,C,D,varargout]=CalcReflCoeff(N,n0,ni)
function [A,B,C,D]=CalcReflCoeff(N,n_out,n_in)
%%
% Function to calculate the reflection coefficient dependent quantities as
% defined in Klose 2005. 
%%
for ii=1:2*N
a=getr(ii,n_out,n_in);                                                
R(ii)=(integral(a,0,pi/2));
end

%Reflection parameters. ref. Klose 
A(1)=-R(1);
B(1)=3*R(2);

% if N>1
    
    A(2)=(-2.25*R(1))+(7.5*R(3))-(6.25*R(5));
B(2)=(15.75*R(2))-(52.5*R(4))+(43.75*R(6));

C(1)=(-1.5*R(1))+(2.5*R(3));
D(1)=(1.5*R(2))-(2.5*R(4));


C(2)=C(1);
D(2)=D(1);
% end

end

function rcoeff=getr(n,n_out,n_in)
%function returns a function handle to evaluate the 'n'th moment of the reflection coefficient for
%given n0 and ni
rcoeff=@refcffni;

    function Rn=refcffni(v1)
    vc=asin(n_out/n_in);
    temp=ones(size(v1));
    %Eqn (5) Klose_2006
    %for defns. of variables chk Klose paper
    ind=find(abs(v1)<vc);
    sinv2=(n_in/n_out)*sin(v1(ind));
    cosv2=sqrt(1-(sinv2).^2);
    temp(ind)=0.5*(((n_in*cosv2-n_out*cos(v1(ind)))./(n_in*cosv2+n_out*cos(v1(ind)))).^2)+...
        0.5*(((n_in*cos(v1(ind))-n_out*cosv2)./(n_in*cos(v1(ind))+n_out*cosv2)).^2);
   Rn=(temp).*(((cos(v1)).^n).*sin(v1));

    end
end