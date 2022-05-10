function [bn,varargout] = get_pfunc_coeff(pfunc,deflag,g,varargin)
% Created on 8th February 2020
%% Obtain Legendre series expansion coefficients for the phase function of
% choice
%%%%%%%%%%% INPUT 
% pfunc - type of phase function. set to 0 - Henyey Greenstein, 
%                                        1 - Modified Henyey Greenstein, 
%                                        2- two term Henyey Greenstein Phase function
%deflag - 0 (default). if 1 - computes the delta-Eddington approximation to
%the phase function 
% g - anisotropy factor
% if pfunc = 1,  varargin{1} = gamma
% if pfunc = 2, varargin{1} = gb , the anisotropy factor for the backscattering
% function (must be negative)
% varargin{2} = alpha , where TTHG = alpha * HG(g) + (1-alpha)*HG(h)
%%%%%%%%%%% OUTPUT
% bn - Legendre series exapnsion coefficients
% if deflag = 1, varragout{1} = f (the delta Eddington parameter)

%%
% Created by Nishigandha Patil
% IIT Kanpur
% nipat@iitk.ac.in, nishi.rpatil@gmail.com
%%
switch(pfunc)
    case 0
        % H-G 
        b(1) = g; b(2) = g^2; b(3) = g^3; b(4) = g^4;
    case 1
        % MHG
        gamma = varargin{1};
        alpha = (gamma-3/5)/(gamma*g-g^2+2/5);
        b(1) = alpha*g; b(2) = alpha*g^2 + (1/5)*(1-alpha); b(3) = alpha*g^3; b(4) = alpha*g^4;
    case 2
        % TTHG
        gb = varargin{2}; alpha = varargin{1};
        b(1) = alpha*g + (1-alpha)*gb; 
        b(2) = alpha*g^2 + (1-alpha)*gb^2; 
        b(3) = alpha*g^3 + (1-alpha)*gb^3; 
        b(4) = alpha*g^4 + (1-alpha)*gb^4; 
end

if (deflag)
    % Delta Eddington approximation 
    % f = b(4); 
    bn(1:3) = (b(1:3)-b(4))/(1-b(4));
    varargout{1} = b(4);
else
    bn(1:3) = b(1:3);
end
                                    
