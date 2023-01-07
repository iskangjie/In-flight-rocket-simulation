function [dis,vel,acc] = Newmark_beta(M,C,K,dof,u,beta,gama,delt_t,Nt,initial)
dis      = zeros(dof,Nt);
vel      = zeros(dof,Nt);
acc      = zeros(dof,Nt);
dis(:,1) = initial(:,1);
vel(:,1) = initial(:,2);
acc(:,1) = M\(u(:,1)-K*dis(:,1)-C*vel(:,1));
alpha_0  = 1/gama/delt_t^2;
alpha_1  = beta/gama/delt_t;
alpha_2  = 1/gama/delt_t;
alpha_3  = 1/2/gama-1;
alpha_4  = beta/gama-1;
alpha_5  = delt_t/2*(alpha_4-1);
alpha_6  = delt_t*(1-beta);
alpha_7  = beta*delt_t;
K_e      = K+alpha_0*M+alpha_1*C;

for k = 1:Nt
    F_e = u(:,k+1)+M*(alpha_0*dis(:,k)+alpha_2*vel(:,k)+alpha_3*acc(:,k))+...
          C*(alpha_1*dis(:,k)+alpha_4*vel(:,k)+alpha_5*acc(:,k));
    dis(:,k+1) = K_e\F_e;
    acc(:,k+1) = alpha_0*(dis(:,k+1)-dis(:,k))-alpha_2*vel(:,k)-alpha_3*acc(:,k);
    vel(:,k+1) = vel(:,k)+alpha_6*acc(:,k)+alpha_7*acc(:,k+1);
end
end