function [bn,varargout] = get_pfunc_coeff(pfunc,deflag,bng,varargin)
% Created on 8th February 2020
% Obtain Legendre series expansion coefficients for the phase function of
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

switch(pfunc)
    case 0
        % H-G 
        b(1) = 1; b(2) = bng; b(3) = bng^2; b(4) = bng^3; b(5) = bng^5;
%         b(1) = bng(1); b(2) = bng(2); b(3) =bng(3); b(4) = bng(4); b(5) = bng(5);
    case 1
        % MHG
        g = varargin{1};
        gamma = varargin{2};
        alpha = (gamma-3/5)/(gamma*g-g^2+2/5);
        b(1) = 1;
        b(2) = alpha*bng; b(3) = alpha*bng^2 + (1/5)*(1-alpha); b(4) = alpha*bng^3; b(5) = alpha*bng^4;
%         b(1) = bng(1);
%         b(2) = alpha*bng(2); b(3) = alpha*bng(3) + (1/5)*(1-alpha); b(4) = alpha*bng(4); b(5) = alpha*bng(5);
    case 2
        % TTHG
        % Not yet corrected for 2D
        g = varargin{1};
        gb = varargin{3}; alpha = varargin{2};
        b(1) = alpha*g + (1-alpha)*gb; 
        b(2) = alpha*g^2 + (1-alpha)*gb^2; 
        b(3) = alpha*g^3 + (1-alpha)*gb^3; 
        b(4) = alpha*g^4 + (1-alpha)*gb^4; 
end

if (deflag)
    % Delta Eddington approximation 
    % f = b(4); 
    bn(1:4) = (b(1:4)-b(5))/(1-b(5));
    varargout{1} = b(5);
else
    bn(1:4) = b(1:4);
end
                                    
