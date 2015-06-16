using Convex
using SCS

# Has helper functions for declaring variables

# Function that runs the offline portion of the risk-averse MPC problem
#
#
#  Note: n_s is the number of states, n_i is the number of inputs and 
#   m is the number of random outcomes per stage.
#
#  Input variable	Size		Description
#  A_scen		n_s,n_s,m	Possible state matrix A values
#  B_scen		n_s,n_i,m	Possible control matrix B values
#  p			1,m		Risk of vertices
#  R_c			n_i,n_i		Cost matrix for constrained control?
#  Q_c			n_s,n_s		Cost matrix for constrained state?
#  R_d			n_i,n_i		Cost matrix for ? control
#  Q_d			n_s,n_s		Cost matrix for ? state
#  T_u			n_i,n_i		Transition matrix for control?
#  T_x			n_s,n_s		Transition matrix for state
#  u_max		1		Control constraint
#  x_max		1		State constraint
#
#  Output variable	size		Description
#  K			
#  Q_out		
#  bar_Q_out		
#  G_out		
#  Y_out		
#  W_N_out		
#  feas_ind		1		Boolean for feasible solution
#  solve_time		1		Time of solution
#
#
function MPC_fn_uncstr(A_scen, B_scen, xi_pts, p, R_c, Q_c,X_init,alpha)

   #Declare variables

   println("Declaring variables")
   feas_ind = 1;
   n_s = size(B_scen[1],1);
   n_i = size(B_scen[1],2);
   m   = size(B_scen,1);
   
  \alpha 
   #### Optimization Variables ####
   
   Q     = fill(Variable(n_s,n_s),m);
   G     = fill(Variable(n_s,n_s),m);
   Y     = fill(Variable(n_i,n_s),m);
   r0    = Variable();

   #Concactenation requires initialization, so we start the loop at 2...
   #Q_diag = hvcat((2,2), Q[1], zeros(n_s,n_s), zeros(n_s,n_s),Q[2]);

   for i = 3:m
      Zi = zeros(n_s*(i-1),n_s);
    #  Q_diag = hvcat((2,2),Q_diag, Zi, Zi',Q[i]);
   end
   Q_diag = eye(n_s*m);

   
   gamma2 = fill(Variable(),m,m,m);
   U2 = fill(Variable(n_i,1),m,m);
   U1 = fill(Variable(n_i,1),m);
   U  = Variable(n_i,1);

   r2 = fill(Variable(),m,m);
   r1 = fill(Variable(),m);


   #### Constants ####
   xtmp = xi_pts';
   #cheat for now.
   Lambda = 0.25*eye(n_s*m,n_s*m);
   Lambda_sqrt = sqrtm(Lambda);
   Sigma  = Lambda;
   #kron(diagm(xtmp[:,1]),eye(n_s));
   Sigma_sqrt = sqrtm(Sigma);
   R_c_sqrt = sqrtm(R_c);
   Q_c_sqrt = sqrtm(Q_c);

   if(size(xi_pts,1) != 1)
      println("Expected xi_pts to be vector. Reconsider Sigma definition")
   end

#cstr = [U ==0];
   #### Constraints ####
   print("Declaring constriants");
#   cstr = [Q[1] > 0.0];
   for i = 2:m
#      cstr+=[Q[i] > 0.0];
   end

   # Enumerate optimality constraints
   for i = 1:m
      b12 = zeros(m*n_s,n_i);
      b13 = zeros(m*n_s,n_s);
      b22 = eye(n_i);
      b23 = zeros(n_i,n_s);
      b33 = eye(n_s);
      b14 = A_scen[1]*G[1] + B_scen[1]*Y[1];
      for i2 = 2:m
         b14 = vcat(b14, A_scen[i2]*G[i2] + B_scen[i2]*Y[i2]);
      end
      b14 = -Sigma_sqrt*b14;

      b24 = -(R_c_sqrt*Y[i]);
      b34 = -(Q_c_sqrt*G[i]);
      b44 = -Q[i] + G[i] + G[i]';

   #    cstr += [hvcat((4,4),Q_diag,b12,b13,b14,
   #                         b12',b22,b23,b24,
#			    b13',b23',b33,b34,
#			    b14',b24',b34',b44)>=0];
#
   end

   # Dynamic constraints

   #Dummy atom for proper typing
   X_tmp = A_scen[1]*X_init + B_scen[1]*U;
   X1 = fill(X_tmp,m);
   X2 = fill(X_tmp,m,m);
   X3 = fill(X_tmp,m,m,m);

   for i1 = 1:m
     X1[i1] = A_scen[i1]*X_init + B_scen[i1]*U;
     for i2 = 1:m
       X2[i1,i2] = A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
       for i3 = 1:m
          X3[i1,i2,i3] = A_scen[i3]*X2[i1,i2] + B_scen[i3]*U2[i1,i2];
          
#	  cstr += [hvcat((2,2), gamma2[i1,i2,i3], X3[i1,i2,i3]',X3[i1,i2,i3], Q[i3]) >= 0];		
       end
     end
   end

   # Risk constraints

#Many "not DCP compliant" complaints
#   sum_exp=0;
#   X_tmp = Variable()+Variable();
#   mean_sum = X_tmp;
#   mean_sum1 = fill(X_tmp, m);
#   mean_sum2 = fill(X_tmp, m,m); 
#   for i1 = 1:m
#      for i2 = 1:m
#        for i3 = 1:m
#           mean_sum2[i1,i2] += p[i3]*max(0,-r2[i1,i2]+gamma2[i1,i2,i3]);
#	end
#	mean_sum2[i1,i2] = r2[i1,i2]+(mean_sum2[i1,i2])/alpha;
#	mean_sum1[i1] += p[i2]*max(0,-r1[i1] + X2[i1,i2]'*Q_c*X2[i1,i2] + U2[i1,i2]'*R_c*U2[i1,i2]+mean_sum2[i1,i2]);
#      end
#      mean_sum1[i1] = r1[i1] + mean_sum1[i1]/alpha;

#      mean_sum += p[i1]*max(0,-r0 + X1[i1]'*Q_c*X1[i1] + U1[i1]'*R_c*U1[i1] + mean_sum1[i1]);
#   end

#   mean_sum = r0 +  mean_sum/alpha;

   gamma1 = X_init'*Q_c*X_init + U'*R_c*U# + mean_sum;
   problem = minimize(gamma1);

   @time solve!(problem,SCSSolver(max_iters=1000000))




end


