function alpha=LineSearch(SetParam,pk,gk,xk,mesh,Src,Det,func,Zold,Md,LS_type,norm_flag)
global SimType
C_armijo=SetParam{1};
beta_armijo=SetParam{2};
tau_armijo=SetParam{3}; %the value by which alpha is decreased. This could be set to a quadrature or cubic rule.
thresh = SetParam{4};
xmin = SetParam{5};
xval_max = SetParam{6};
alpha_min = 1e-3;
alpha=1;
stop=0;
m=squeeze(pk)'*gk;
t_armijo=C_armijo*m;

if LS_type ==1
    % Backtracking line search
    while stop==0
        xnew = xk + alpha*pk;
        xnew(xnew<=thresh*max(xnew)) = xmin;
        xnew(xnew>xval_max)=xval_max;
%         Param{SimType+1}=xnew;
        mesh.opt.muaxf =xnew;
        [Rnew,~,~] = eval_residual(xnew,Md,func,mesh,Src,Det,norm_flag);
        Znew=0.5*(Rnew'*Rnew);
        if (Znew<=(Zold+alpha*beta_armijo*t_armijo))
            stop=1;
        else if alpha<=alpha_min;
                stop=1;
                alpha=0;
            else
                alpha=tau_armijo*alpha;
            end
        end
    end
else
    % Quadratic-Cubic Line Search
    
    xnew = xk + alpha*pk;
    xnew(xnew<thresh) = xmin;
%     xnew(xnew>xval_max)=xval_max;
%     Param{SimType+1}=xnew;
    mesh.opt.muaxf = xnew;
    [Rnew,~,~] = eval_residual(xnew,Md,func,mesh,Src,Det,norm_flag);
    Znew=0.5*(Rnew'*Rnew);
    if (Znew<=(Zold+alpha*beta_armijo*t_armijo))
        stop=1;
    else
        % Try quadratic fit
        alpha_q = -(norm(gk)*alpha^2)/(2*(Znew-Zold-alpha*norm(gk)));
        if (abs(alpha_q-alpha)<=1e-2)
            alpha_q = alpha/2; % Reinitialise alpha
        end
        xnew = xk + alpha_q*pk;
        xnew(xnew<thresh) = xmin;
%         xnew(xnew>xval_max)=xval_max;
%         Param{SimType+1}=xnew;
        mesh.opt.muaxf = xnew;
        Rq = eval_residual(xnew,Md,func,mesh,Src,Det,norm_flag);
        Zq=0.5*(Rq'*Rq);
        if (Zq<=(Zold+alpha*beta_armijo*t_armijo))
            stop=1;
        else
            loop =1;
            while(loop)
                % Try cubic fit till you find appropriate alpha
                div = (1/((alpha_q^2*alpha^2)*(alpha_q-alpha)));
                a= ((Zq - Zold -norm(gk)*alpha_q)*alpha^2 - (Znew-Zold-norm(gk)*alpha)*alpha_q^2)*div;
                b = (-(Zq - Zold -norm(gk)*alpha_q)*alpha^3 + (Znew-Zold-norm(gk)*alpha)*alpha_q^3)*div;
                alpha_c=(-b+sqrt(b^2 - 3*a*norm(gk)))/(3*a);
                if (abs(alpha_c-alpha_q)<=1e-2)
                    alpha_c=alpha_q/2;
                    if abs(alpha_c)<=1e-3
                        
                        alpha=0;
                    end
                end
                xnew = xk + alpha_c*pk;
                xnew(xnew<thresh) = xmin;
%                 xnew(xnew>xval_max)=xval_max;
%                 Param{SimType+1}=xnew;
                mesh.opt.muaxf = xnew;
                Rc = eval_residual(xnew,Md,func,mesh,Src,Det,norm_flag);
                Zc=0.5*(Rc'*Rc);
                if (Zc<=(Zold+alpha*beta_armijo*t_armijo))
                    loop = 0;
                    alpha = alpha_c;
                else
                    Znew=Zq;
                    Zq=Zc;
                    if (abs(alpha)>1)||~isreal(alpha_q)
                        loop = 0;
                        alpha = 0;
                    else
                        alpha=alpha_q;
                        alpha_q=alpha_c;
                    end
                end
            end
        end 
    end
end
