function pk = folm_step(Jacobian,Scale_mat,Residual,lambda,L)
[Qj,Rj,Pj] = qr([Jacobian*Scale_mat;sqrt(lambda)*L]);
pk=Scale_mat*Pj*(-Rj\Qj'*[Residual;zeros(size(L,1),1)]);
end