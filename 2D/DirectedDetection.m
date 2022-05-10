function Jd=DirectedDetection(N,nfibre,ni,AOD)
%% Calculate J -type coefficients for partial current measurement when collection is
% in a pre-specified direction. 
%%%%%%%%%%%% INPUT
% N - order of SPN approximation
% n0 - refractive index of outside medium / fiber
% ni - refractive index of the medium
% AOD - angle of detection
%%%%%%%%%%%% OUTPUT
% J - [(N+1) X 1] vector containing J-type coefficients. 
% Refer Klose et al for details. 

%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
AOD=AOD*pi/180;
vc=asin(nfibre/ni); %Critical angle of incidence from inside the medium
v1=asin((nfibre/ni)*sin(AOD));
if v1<vc
    sinv2=sin(AOD);
    cosv2=sqrt(1-(sinv2).^2);
    cosv1=cos(v1);
R=0.5*(((ni*cosv2-nfibre*cosv1)./(ni*cosv2+nfibre*cosv1)).^2)+...
        0.5*(((ni*cosv1-nfibre*cosv2)./(ni*cosv1+nfibre*cosv2)).^2);
else
    R=1;
end
%Define Legendre polynomials 
P0=1;
P1=cos(AOD);
P2=0.5*(3*P1^2-1);
P3=0.5*(5*P1^3-3*P1);

Jd(1)=(1-R)*P1*P0/2;
Jd(3)=(1-R)*P1*(-P0/3+P2*5/6);
Jd(2)=(1-R)*P1*P1;
Jd(4)=(1-R)*P1*P3;

if N>3
    P4=0.25*((35/2)*P1^4-15*P1^2+3/2);
    P5=0.25*((63/2)*P1^5-35*P1^3+(15/2)*P1);

    Jd(5)=(1-R)*P1*(4*P0/15-2*P2/3+9*P4/10);
    Jd(6)=(1-R)*P1*P5;
    if N>5
        P6=(1/16)*(231*P1^6-315*P1^4+105*P1^2-5);
        P7=(1/16)*(429*P1^7-693*P1^5+315*P1^3-35);
        Jd(7)=(1-R)*P1*(-8*P0/35+8*P2/14-27*P4/35+13*P6/14);
        Jd(8)=(1-R)*P1*P7;
    end
end
