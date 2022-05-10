function [lambda,varargout] = update_lambda(k,lambda,rho,rho_accept,rho_dec,update_type,varargin)

switch (update_type)
    case 1
          % Continuous updation strategy from Madsen
          q=varargin{1};
        dec_factor = 1/3;
        inc_factor = 2;
        if (rho>=rho_accept)
            lambda(k+1) = lambda(k)*min(dec_factor,1-(dec_factor-1)*(1-2*rho)^3);
            q = inc_factor;
        elseif rho<rho_accept
            lambda(k) = q*lambda(k);
            q=inc_factor*q;
        end
        varargout{1}=q;
    case 2
        % Basic updation strategy 2:3
        dec_factor = 1/3;
        inc_factor = 2;
        if (rho>=rho_accept)&&(rho>=rho_dec)
            lambda(k+1) = dec_factor*lambda(k);
        elseif (rho>=rho_accept)&&(rho<rho_dec)
            lambda(k+1) = lambda(k);
        else
            lambda(k)=inc_factor*lambda(k);
        end
        
    case 3
        % Basic updation strategy 3:5
         dec_factor = 1/5;
        inc_factor = 3;
        if (rho>=rho_accept)&&(rho>=rho_dec)
            lambda(k+1) = dec_factor*lambda(k);
        elseif (rho>=rho_accept)&&(rho<rho_dec)
            lambda(k+1) = lambda(k);
        else
            lambda(k)=inc_factor*lambda(k);
        end
        
        

end