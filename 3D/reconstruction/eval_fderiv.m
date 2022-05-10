function Jacobian = eval_fderiv(x,func_deriv,mesh,Src,Det,Jd,FwData,norm_data,options)
%% Created on 27 December 2019
global Dimension N AssemblyFile SimType fcs_flag JacType pfuncflag

xk(1,1,:)=x;
if (SimType) 
    mesh.opt.muaxf=xk;
else
    mesh.opt.muaxi = xk;
end
NofSources = size(Src,2);

switch(norm_data.flag)
    case 0
        % No normalisation
        Jac = func_deriv(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,options,FwData);
        for s=1:NofSources
            F(:,s)=Jd{s}(:,SimType+1);
        end
        F = F(:);
        % Uncomment for log scaling
        if isreal(F)
            Jacobian=(diag(F)\Jac);%imag(diag(F)\J)];
        else
            Jacobian=[real(diag(F)\Jac);imag(diag(F)\Jac)];
        end
%          if isreal(F)
%             Jacobian=Jac;
%         else
%             Jacobian=[real(Jac);imag(Jac)];
%         end
        
    case 1
        % Normalisation with excitation data. Only for SimType =1;
        [JacX,JacM] = func_deriv(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,options,FwData);
        for s=1:NofSources
            FX(:,s)=Jd{s}(:,SimType);
            FM(:,s)=Jd{s}(:,SimType+1);
        end
        FX = FX(:);
        FM=FM(:);
        if isreal(FX)
            Jacobian=(-diag(FX)\JacX) + (diag(FM)\JacM);%imag(diag(F)\J)];
        else
            Jacobian=[real(-diag(FX)\JacX + diag(FM)\JacM);imag(-diag(FX)\JacX + diag(FM)\JacM)];
        end
        
    case 2
        % Normalisation with respect to ref detector
        Jac = func_deriv(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,options,FwData);
%         if norm_data.bscat
            % data in backscattering mode
           
            ndet = norm_data.ndet;
            Jacobian = [];
            J_real=[]; J_imag=[];
            for s=1:NofSources
                 ref_det = norm_data.ref_det(s);
                F=Jd{s}(:,SimType+1);
                Fref = Jd{s}(ref_det,SimType+1);
                if isreal(F)
                    Jacobian=[Jacobian;(diag(F)\Jac((s-1)*ndet+1:s*ndet,:)) - repmat(Fref\Jac((s-1)*ndet + ref_det,:),[ndet,1])];%imag(diag(F)\J)];
                else
                    temp = (diag(F)\Jac((s-1)*ndet+1:s*ndet,:)) - repmat(Fref\Jac((s-1)*ndet + ref_det,:),[ndet,1]);
                    J_real=[J_real; real(temp)];
                    J_imag =[J_imag; imag(temp)];
                end
            end
            if ~isreal(F)
            Jacobian = [J_real;J_imag];
            end
%         else
%             % data on all sides
%             ndet = norm_data.ndet;
%             Jacobian = [];
%             J_real=[]; J_imag=[];
%             for s=1:NofSources
%                 for side =1:4
%                     ref_det = (side-1)*ndet + norm_data.ref_det;
%                     F=Jd{s}((side-1)*ndet+1:side*ndet,SimType+1);
%                     Fref = Jd{s}(ref_det,SimType+1);
%                     jac_ind = (s-1)*4*ndet + (side-1)*ndet +1: (s-1)*4*ndet + side*ndet;
%                     if isreal(F)
%                         Jacobian=[Jacobian;(diag(F)\Jac(jac_ind,:)) - repmat(Fref\Jac((s-1)*4*ndet + ref_det,:),[ndet,1])];%imag(diag(F)\J)];
%                     else
%                         temp = (diag(F)\Jac(jac_ind,:)) - repmat(Fref\Jac((s-1)*4*ndet + ref_det,:),[ndet,1]);
%                         J_real=[J_real; real(temp)];
%                         J_imag =[J_imag; imag(temp)];
%                     end
%                 end
%             end
%             Jacobian = [J_real;J_imag];
%         end
        
        
end

% if ~fcs_flag
%     Jac = func_deriv(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% else
%     Jac = func_deriv(Dimension,N,SimType,JacType,mesh,Src,Det,AssemblyFile,pfuncflag,FwData);
% end
% for s=1:NofSources
%             F(:,s)=Jd{s}(:,SimType+1);
% end
% F = F(:);
% if isreal(F)
%     Jacobian=(diag(F)\Jac);%imag(diag(F)\J)];
% else
%     Jacobian=[real(diag(F)\Jac);imag(diag(F)\Jac)];
% end