function [d,f] = d_f_calculation(L,A,B,C,u_max,u_min,M,Ap,Bp,Cp,te,d_on,f_on)

%%%%%%%Å@dÅCÉ’ÇãÅÇﬂÇÈ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(d_on == 1)
    psi0 = zeros(L+1,1);
    for K = 0 : L
        psi0(K+1,1) = abs(C*(A+B*C)^K*B);
    end
    psi = sum(psi0);
    psi_max_L =  psi;
    psi_min_L = -psi;

    [T, lam] = eig(A+B*C);
   
    lambda = max(max(abs(lam)));
    
    psi_hat = (abs(C*T)*abs(T\B)*abs(lambda)^L) / (1-abs(lambda));

    d = (u_max-u_min) / (M-psi_hat-(psi_max_L-psi_min_L)/2);
end
%%%%%%%f=E(Q)ÇãÅÇﬂÇÈ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(f_on == 1)
    F = zeros(te,1); 
    A_bar = [Ap                Bp*C;
             zeros(size(B*Cp)) A+B*C];
    B_bar = [Bp;
             B];
    C_bar = [Cp zeros(size(C))];   
   
    for t = 1 : te
        F(t,1) = abs(C_bar*A_bar^t*B_bar);
    end
    
    f = norm((abs(Cp*Bp)+sum(F)),Inf)*(d/2);

    h1 = max(abs(eig(A+B*C))) - 1;
    if(d<=0||h1>=0)
        f1 = 10+1*rand;
        if(f1 > f)
            f = f1;
        end
    end
    
else
    f=1000000;
end

end