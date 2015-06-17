
include("MPC_fn_offline.jl")
function MPC_test()

   #Declare necessary variables
   n_rand = 4;
   n_state = 2;
   n_input = 2;

   A_scen = fill(eye(n_state),n_rand);
   B_scen = fill(eye(n_state),n_rand);

   xi_pts = 0.25*ones(80,n_rand);
   p      = 1/n_rand*ones(1,n_rand);
   R_c    = 0.0001*eye(n_input);
   R_d    = 0.0001*eye(n_input);
   Q_c    = 0.25*eye(n_state);
   Q_d    = 0.0001*eye(n_state);
   T_u    = 0.1*eye(n_input);
   T_x    = 0.1*eye(n_state);
   u_max  = 50;
   x_max  = 50;


   X_init = ones(n_state,1);
   alpha = 0.5


   println("Running nominal test");

   MPC_fn_uncstr(A_scen, B_scen, xi_pts, p, R_c, Q_c, X_init, alpha);



end
