function pk = folm_step_gmres(Jacobian,Scale_mat,gk,lambda,L)

gfunc = @get_hx;

pk = gmres(gfunc,gk);

    function Hx =  get_hx(x)
    Hx = Jacobian*Scale_mat*x;
    Hx = (Jacobian*Scale_mat)'*Hx + lambda*(L'*L)*x;
    end

end


