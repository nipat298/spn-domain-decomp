function [residual,Jd,FwData] = eval_residual(x,Md,func,mesh,Src,Det,norm_data,options)
%% 27 December 2019
% Modified on 19 September 2020 - supports normalisation with reference detector
% as well as with excitation measurements
% include normalisation type here in future
global Dimension N AssemblyFile SimType fcs_flag evalFluence probType pfuncflag

xk(1,1,:)=x;
if (SimType)
    mesh.opt.muaxf=xk;
else
    mesh.opt.muaxi = xk;
end
NofSources = size(Src,2);

if ~fcs_flag
    [Jd,FwData]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,options);
else
    [~,Jd,FwData]=func(Dimension,N,SimType,mesh,Src,Det,AssemblyFile,options);
end

if (norm_data.flag==1)&&(SimType==1)
    for s=1:NofSources
        Jd_n=Jd{s}(:,SimType+1)./Jd{s}(:,SimType);
        if isreal(Jd_n)
            residual(:,s)=log(Jd_n)-Md{s};
        else
            residual_real(:,s)=real(log(Jd_n))-Md.real{s};
            residual_imag(:,s)=imag(log(Jd_n))-Md.imag{s};
        end
    end
elseif norm_data.flag==2
    %     if norm_data.bscat
    % backscatter mode - data only on one side
    for s=1:NofSources
        Jd_n=Jd{s}(:,SimType+1)/Jd{s}(norm_data.ref_det(s),SimType+1);
        if isreal(Jd_n)
            residual(:,s)=log(Jd_n)-Md{s};
        else
            residual_real(:,s)=real(log(Jd_n))-Md.real{s};
            residual_imag(:,s)=imag(log(Jd_n))-Md.imag{s};
        end
    end
    %     else
    %         % det on all 4 sides. need to normalise for each side separately
    %         ndet = norm_data.ndet;
    %
    %         for s=1:NofSources
    %             for side=1:4
    %                 ref_det = (side-1)*ndet + norm_data.ref_det;
    %                 det_ind = (side-1)*ndet+1:side*ndet;
    %                 Jd_n=Jd{s}(det_ind,SimType+1)/Jd{s}(ref_det,SimType+1);
    %                 if isreal(Jd_n)
    %                     residual(det_ind,s)=log(Jd_n)-Md{s}(det_ind);
    %                 else
    %                     residual_real(det_ind,s)=real(log(Jd_n))-Md.real{s}(det_ind);
    %                     residual_imag(det_ind,s)=imag(log(Jd_n))-Md.imag{s}(det_ind);
    %                 end
    %             end
    %         end
    %
    %     end
else
    % No normalisation
    for s=1:NofSources
        Jd_n=Jd{s}(:,SimType+1);
                if isreal(Jd_n)
                    residual(:,s)=log(Jd_n)-Md{s};
                else
                    residual_real(:,s)=real(log(Jd_n))-Md.real{s};
                    residual_imag(:,s)=imag(log(Jd_n))-Md.imag{s};
                end
%         if isreal(Jd_n)
%             residual(:,s)=(Jd_n)-Md{s};
%         else
%             residual_real(:,s)=real((Jd_n))-Md.real{s};
%             residual_imag(:,s)=imag((Jd_n))-Md.imag{s};
%         end
    end
end
if isreal(Jd_n)
    residual = residual(:);
else
    residual = [residual_real(:);residual_imag(:)];
end
