function ReconData = ResidualMinAlgo_optim(x0,Md,func,func_deriv,mesh,Src,Det,nls_method,LS_type,update_type,L,norm_data,solver_options,varargin)
global thresh_x xval_min evalFluence probType SimType
%%%%%%%%%% INPUT
% xk - input parameter
% Md - measurements
% func - func call to fwdsolver
% FwdCallData - supplementary data need in call to func
% nls_method = method for solving the non-linear least squares problem
%               '1' - First order Levenberg Marquardt
%               '2' - Second order Levenberg Marquardt
%               '3' - Levenberg-marquardt with Geodesic acceleration
%               '4' - Tensor Newton with first order inner iteration ref. Gould, RALfit documentation
% LS_type - type of globalisation strategy
%               '1' - backtracking line search
%               '2' - quad-cube line search
%               '3' - none
%update_type - Choice of update strategy
%               '1' - Madsen
%               '2' - simple 2:3
%               '3' - simple 3:5
%%
%% Max counts and tolerances
MAXITER= 50; % number of maximum iterations
max_count_static_cost = 3; % number of counts for which ther is no substantial change in cost function

tol_gradient = 1e-2; % min allowed gradient value
tol_step = 1e-2;    % min step size allowed
tol_cost = 1e-4;    % minimum cost function val

rho_accept = 0.01;
rho_dec = 0.3;
thresh_x = 0;
xval_min =0.031;
xval_max = 1;

support_constraint = 1;
if (support_constraint)
    support_ymin = -0.8/2;
    support_ymax = +0.8/2;
    support_xmin=0.5-0.8/2;
    support_xmax=0.5+0.8/2;
    
     pF1= mesh.nodes((mesh.nodes(:,1)>=support_xmin)&(mesh.nodes(:,1)<=support_xmax)...
            &(mesh.nodes(:,2)>=support_ymin)&(mesh.nodes(:,2)<=support_ymax),:);
        for i=1:size(pF1,1)
            Nds_F1(i)=find((mesh.nodes(:,1)==pF1(i,1))&(mesh.nodes(:,2)==pF1(i,2)));
        end
        BN=ismember(mesh.tri,Nds_F1);
        support_mask=diag((sum(BN,2)==3));
else
    support_mask = 1; 
end
%% FLAGS
stop = 0;
scale =1;
if (LS_type==1)||(LS_type ==2)
    LS = 1;
else
    LS = 0;
end
% reg_update_flag = 0; % Set to 1 to update regularisation parameter

%%
Zk = zeros(MAXITER,1);
lambda = zeros(MAXITER,1);
xk = zeros(length(x0),MAXITER);
rho = zeros(MAXITER,1);

q=2;
k=1;
count = 0;
lambda(1) = 10;
eval=1;
%%
xk(:,k)=x0;
while stop==0
    if eval ==1
        fprintf('************************** Iteration number %d  *******************\n',k);
        [Rk,Jd,FwData] = eval_residual(xk(:,k),Md,func,mesh,Src,Det,norm_data,solver_options);
        Zk(k) = 0.5*(Rk'*Rk);
        
        if Zk(k)<10*tol_cost
            stop =4;
            disp('Exiting on convergence...');
        end
        
%         if ~norm_flag
            Jk = eval_fderiv(xk(:,k),func_deriv,mesh,Src,Det,Jd,FwData,norm_data,solver_options);
%         else
%             Jk = eval_fderiv_mx(xk(:,k),func_deriv,mesh,Src,Det,Jd,FwData);
%         end
        if scale ==1
            Scale_mat = diag(1./(sqrt(sum(Jk.*Jk,1))));
        else
            Scale_mat = diag(ones(size(Jk,2),1));
        end
        gk = Scale_mat'*support_mask'*Jk'*Rk;
        
%         if nls_method~=1
%             if isempty(varargin)
%                 freeze_at = MAXITER;
%             else
%                 freeze_at = varargin{1};
%             end
%             if k <= freeze_at
%                 clear Hess
%                 Hess = eval_sderiv_optim(xk(:,k),func,Param);
%                 beep
%             end
%         end
    end
    if norm(gk,'inf')<=tol_gradient
        stop=1;
        %Returning because slope is too small
        disp('Gradient is too small')
    else
        
        switch (nls_method)
            case 1
                pk = folm_step(Jk*support_mask,Scale_mat,Rk,lambda(k),L*support_mask);
                
%                 pk = folm_step_gmres(Jk,Scale_mat,-gk,lambda(k),L);
%             case 2
%                 pk = solm_step_optim(Jk,Hess,Scale_mat,Rk,lambda(k),L);
%             case 3
%                 pk = galm_step_optim(Jk,Hess,Scale_mat,Rk,lambda(k),L);
%             case 4
%                 
%                 %                  [pk,mk] = tnlm_step_FO_optim1_1(xk(:,k),Jk,Hess,Rk,lambda(k),L,update_type);
%                 [pk,mk] = tnlm_step_FO_optim2_2(xk(:,k),Jk,Hess,Rk,lambda(k),L,update_type);
%                 %                 [pk,mk] = tnlm_step_FO_optim3(xk(:,k),func,Md,Param,Jk,Hess,Rk,lambda(k),L,update_type);
                %                 [pk,mk] = tnlm_step_SO_optim2(xk(:,k),Jk,Hess,Rk,lambda(k),L,update_type);
                
        end
        if norm(squeeze(pk))<=tol_step*(norm(squeeze(xk(:,k)))+tol_step)
            if norm(squeeze(pk))~=0
                stop=2;
                %Too small change
                disp('Step size too small')
            else
                lambda = update_lambda(k,lambda,0,rho_accept,rho_dec,update_type);
            end
        else
            if LS ==1
                LS_Param{1}=0.1; %sets the value for C_armijo
                LS_Param{2}=1; %sets the value for beta_armijo
                LS_Param{3}=0.5;  %sets the value for tau_armijo
                LS_Param{4} = thresh_x; % Threshold value for xk. Set to Zero if no thresholding is desired
                LS_Param{5} = xval_min; % Minimum value of xk to set after thresholding.
                LS_Param{6} = xval_max;
                alpha=LineSearch(LS_Param,pk,gk,xk(:,k),mesh,Src,Det,func,Zk(k),Md,LS_type,norm_data);
            else
                alpha =1;
            end
            if alpha>0
                x_new=xk(:,k)+alpha*pk;
                x_new(x_new<=thresh_x*max(x_new))= xval_min;%min(xval_min,thresh_x*max(x_new));
                x_new(x_new>xval_max)=xval_max;
                
%                if (support_constraint)
%                     pF1= mesh.nodes((mesh.nodes(:,1)>=support_xmin)&(mesh.nodes(:,1)<=support_xmax)&(mesh.nodes(:,2)>=support_ymin)&(mesh.nodes(:,2)<=support_ymax),:);
%                      for i=1:size(pF1,1)
%                             Nds_F1(i)=find((mesh.nodes(:,1)==pF1(i,1))&(mesh.nodes(:,2)==pF1(i,2)));
%                      end
%                     BN=ismember(mesh.tri,Nds_F1);
%                     tF1=(sum(BN,2)==3);
%                     x_new(~tF1) = xval_min;
%                end
%                 x_new(x_new<xval_min)=xval_min;
                %         x_new = xk+pk;
                [Rnew,~,~] = eval_residual(x_new,Md,func,mesh,Src,Det,norm_data,solver_options);
                Znew = 0.5*(Rnew'*Rnew);
                if nls_method==4
                    rho(k) = (Zk(k)-Znew)/abs(Zk(k)-mk);
                else
                    rho(k) = (Zk(k)-Znew)/(0.5*squeeze(pk)'*(lambda(k)*(L'*L)*squeeze(pk)-gk));
                end
                fprintf('lambda = %f \n rho = %f \n Cost function val = %f \n', lambda(k) ,rho(k), Znew);
            else
                
                rho(k) =0;
            end
            
%             beep
%             pause(2)
%             beep
            if update_type ==1
                [lambda,q] = update_lambda(k,lambda,rho(k),rho_accept,rho_dec,update_type,q);
            else
                lambda = update_lambda(k,lambda,rho(k),rho_accept,rho_dec,update_type);
            end
            
            if rho(k)>rho_accept
                k=k+1;
                xk(:,k)=x_new;
                eval=1;
            else
                eval=0;
            end
            
            if k>1
                if abs(Zk(k-1)-Znew)<tol_cost
                    count = count+1;
                    if count ==max_count_static_cost
                        stop = 3;
                        disp('Cost function unchanged for 5 iterates....Consider restarting the algo')
                    end
                else
                    count =0;
                end
            end
            
            if k>MAXITER
                stop=5;
                disp('Exiting on MAXITER')
            end
            
            if lambda(k)>1e4
                stop = 2;
                disp('Lambda too large')
            end
            
            fprintf('norm(gk) = %f, norm(rk) = %f, Proj_error = %f \n',norm(gk,'inf'),norm(Rk),norm(gk,'inf')/norm(Rk))
            if (norm(gk,'inf')/norm(Rk))<tol_cost
                stop = 5;
                disp('Exiting on convergence...');
            end        
        end
    end
end
ReconData{1}=xk;
ReconData{2}=Zk;
ReconData{3}=lambda;
ReconData{4}=rho;

end





