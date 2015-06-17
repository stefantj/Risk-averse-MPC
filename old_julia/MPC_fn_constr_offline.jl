using Convex
using SCS

# Has helper functions for declaring variables
include("MPC_helpers.jl");



# Function that runs the offline portion of the risk-averse MPC problem
#
#
#  Note: n_s is the number of states, n_i is the number of inputs and 
#   m is the number of random outcomes per stage.
#
#  Input variable	Size		Description
#  A_scen		n_s,n_s,m	Possible state matrix A values
#  B_scen		n_s,n_i,m	Possible control matrix B values
#  xi_pts		1,m		Vertices of risk polytope
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
function MPC_fn_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c, R_d, Q_d, T_u, T_x, u_max, x_max)


   #Declare variables

   println("Declaring variables")
   feas_ind = 1;
   n_s = size(B_scen,1);
   n_i = size(B_scen,2);
   m   = size(B_scen,3);
   

   #### Optimization Variables ####
   
   #Not really sure what W_N is, but it looks like a slack variable
   W_N = Variable(n_s,n_s);

   # n-d variable in julia
   Q     = fill(Variable(n_s,n_s),m);
   bar_Q = fill(Variable(n_s,n_s),m);
   G     = fill(Variable(n_s,n_s),m);
   Y     = fill(Variable(n_i,n_s),m);
   
   #Concactenation requires initialization, so we start the loop at 2...
   Q_diag = Q[1];
   bar_Q_diag = bar_Q[1];
   A_all = A_scen[:,:,1];
   B_all = B_scen[:,:,1];

   for i = 2:m
      # Reshape A_scen, B_scen
      A_all = vcat(A_all, A_scen[:,:,i]);
      B_all = vcat(B_all, B_scen[:,:,i]);

      Zi = zeros(n_s*(i-1),n_s);
      Q_diag = hvcat((2,2),Q_diag, Zi, Zi',Q[i]);
      bar_Q_diag = hvcat((2,2),bar_Q_diag,Zi,Zi',bar_Q[i]);
   end

   #### Constants ####
   Lambda = kron(diagm(p[:]),eye(n_s)); #check this.
   Lambda_sqrt = sqrtm(Lambda);
   Sigma  = kron(diagm(xi_pts[:]),eye(n_s));
   Sigma_sqrt = sqrtm(Sigma);
   R_d_sqrt = sqrtm(R_d);
   R_c_sqrt = sqrtm(R_c);
   Q_d_sqrt = sqrtm(Q_d);
   Q_c_sqrt = sqrtm(Q_c);

   if(size(xi_pts,1) != 1)
      println("Expected xi_pts to be vector. Reconsider Sigma definition")
   end

   #predeclare some helpful matrices
   Z_ns = zeros(n_s,n_s);
   Z_ni = zeros(n_i,n_i);
   Z_mns_ni = zeros(m*n_s,n_i);
   Z_mns_ns = zeros(m*n_s,n_s);
   Z_ni_ns = zeros(n_i,n_s);
   I_ns = eye(n_s);
   I_ni = eye(n_i);

   #### Constraints ####
   print("Declaring constriants");
   cstr = [W_N > 0];
   cstr += [Q_diag > 0];
   cstr += [bar_Q_diag > 0];

   # Enumerate optimality constraints
   for i = 1:m
      b12 = zeros(m*n_s,n_i);
      b13 = zeros(m*n_s,n_s);
      b22 = eye(n_i);
      b23 = zeros(n_i,n_s);
      b33 = eye(n_s);

      b14 = -Sigma_sqrt*(A_all*G[i] + B_all*Y[i]);
      b24 = -R_c_sqrt*Y[i];
      b34 = -Q_c_sqrt*G[i];
      b44 = -Q[i] + G[i] + G[i]';

      cstr += [hvcat((4,4),Q_diag,b12,b13,b14,
                         b12',b22,b23,b24,
			 b13',b23',b33,b34,
			 b14',b24',b34',b44) >= 0];

#      b14 = -Lambda_sqrt*(A_all*G[i] + B_all*Y[i]);
#      b24 = -R_d_sqrt*Y[i];
#      b34 = -Q_d_sqrt*G[i];
#      b44 = -bar_Q[i] + G[i] + G[i]';
#
#      cstr += [hvcat((4,4),bar_Q_diag,b12,b13,b14,
#                                 b12',b22,b23,b24,
#			         b13',b23',b33,b34,
#			         b14',b24',b34',b44) >= 0];
#

#      b11 = G[i]+G[i]' - W_N;
#      cstr += [hvcat((2,2),b11, (T_u*Y[i])',
#                      T_u*Y[i], u_max^2*eye(n_i)) >= 0];
#      for j = 1:m
#         b12 = (A_scen[:,:,j]*G[i] + B_scen[:,:,j]*Y[i])';
#         cstr += [hvcat((2,2),b11,b12, b12', W_N) >= 0];
#	  cstr += [hvcat((2,2),b11,b12*(T_x'), T_x*(b12'), x_max^2*eye(n_s)) >= 0];
#
#      end
   end


   objective = -logdet(Q_diag)-logdet(bar_Q_diag)-logdet(W_N);
   problem = minimize(objective, cstr);

   @time solve!(problem,SCSSolver(max_iters=1000000))




end


