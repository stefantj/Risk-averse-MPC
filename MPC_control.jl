using Convex
using SCS
using Debug

# This file has functions for executing risk-averse MPC efficiently
#
# MPC_test	Tests control problem with the MPC functions
# 
# MPC_look_2	Look-ahead of 2
# MPC_look_3	Look-ahead of 3
# MPC_look_4	Look-ahead of 4
# MPC_look_5	Look-ahead of 5
# 
# MPC_offline	Pre-solves terminal constraints as suggested in paper
#
# 


# Function that computs Q offline
function MPC_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs

  # Decision variables for the LMI that describes optimality
  Q = [Semidefinite(n_s) for i=1:m];
  G = [Variable(n_s,n_s) for i=1:m];
  Y = [Variable(n_i,n_s) for i=1:m];

  # Useful constants that are used a lot
  R_inv    = inv(R_c);
  R_c_sqrt   = sqrtm(R_c);
  Q_c_sqrt   = sqrtm(Q_c);

  # Problem constraints

  # Q > 0
  cstr = Q[1] > 0.0000001*eye(n_s);
  for i = 2:m
     cstr+= Q[i] > 0.000000001*eye(n_s);
  end

  # Optimality constraint
  #  note: using the LMI from the published paper versus the sample code
  for i = 1:m
      Q_diag = hvcat((2,2),Q[i],zeros(n_s,n_s),zeros(n_s,n_s),Q[i]);
      for j = 3:m
        Zi = zeros(n_s*(j-1), n_s);
        Q_diag = hvcat((2,2),Q_diag,Zi,Zi',Q[i]);
      end

      b12 = zeros(m*n_s,n_i);
      b13 = zeros(m*n_s,n_s);
      b22 = R_inv;
      b23 = zeros(n_i,n_s);
      b24 = -(Y[i])
      b33 = eye(n_s);
      
      #sigma matrix not necessary, since for our example everything is constant
      b14 = -0.5*(A_scen[1]*G[1] + B_scen[1]*Y[1]);
      for j = 2:m
        b14 = vcat(b14, -0.5*(A_scen[j]*G[j] + B_scen[j]*Y[j]));
      end

      b34 = -Q_c_sqrt*G[i];
      b44 = -Q[i] + G[i] + G[i]';

      cstr+= [hvcat((4,4,4,4), Q_diag, b12, b13, b14,
                           b12', b22, b23, b24,
			   b13', b23', b33, b34,
			   b14', b24', b34', b44) >= 0];
  end

  # Time to solve the problem!
  
  objective = -logdet(Q[1]);
  for i = 2:m
    objective -= logdet(Q[i]);
  end

  problem = minimize(objective, cstr);
  solve!(problem,SCSSolver(verbose=false));
  return Q, problem
end

function MPC_look_2(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init, Q)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs

  # Control inputs
  U  = Variable(n_i,1);
  U1 = [Variable(n_i,1) for i=1:m];

  # Trajectories
  X1 = [Variable(n_s,1) for i=1:m];
  X2 = [Variable(n_s,1) for i=1:m, j = 1:m];

  # Epigraph variables
  gamma1 = Variable();
  gamma2 = [Variable() for i=1:m,j=1:m];

  # Problem constraints
  cstr = []
  # System dynamics
  for i1 = 1:m
    cstr += X1[i1] == A_scen[i1]*X_init + B_scen[i1]*U;
    for i2 = 1:m
      cstr += X2[i1,i2] == A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
      # Final action constraint (a byproduct of using epigraph to minimize)
      cstr += hvcat((2,2),gamma2[i1,i2],X2[i1,i2]',X2[i1,i2],Q[i2]) >= 0
    end
  end

  # Cost epigraph constraints 
  expected_cost = quad_form(U,R_c);
  
  # Move these out at some point.
  r0 = 0.25;
  max_risk = 0.25;
  alpha = 0.25;

  r1 = [Variable() for i=1:m]

  MS = 0;  

  for i1 = 1:m
    # Cost of first action
    c_ij = quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c);
    MS2 = 0;
    for i2 = 1:m
      c_ij += gamma2[i1,i2];
      expected_cost += c_ij*p[i1]*p[i2];

      # Condition 2 - check Q??
      cstr += hvcat((2,2),gamma2[i1,i2], X2[i1,i2]',X2[i1,i2],Q[i2]) >= 0;
      MS2 = p[i2]*max(0,-r1[i1] + gamma2[i1,i2]);
    end
    MS2 = r1[i1] + MS2/alpha;
    MS = p[i1]*max(0,-r0 + quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c) + MS2);
  end
  MS = r0 + MS/alpha;

  cstr += gamma1 == quad_form(U,R_c) + MS; #Note that this does not include the initial cost X_init'*Q_C*X_init
  cstr += expected_cost <= max_risk;

  # Time to solve the problem!
  problem = minimize(gamma1,cstr);
  solve!(problem, SCSSolver(verbose=false))

  # TODO: second argument should be expected cost. 'expected_cost' is addition_atom, hence has no value.
  return U.value, 0.0, problem 
end

function MPC_look_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init, Q)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs

  # Control inputs
  U  = Variable(n_i,1);
  U1 = [Variable(n_i,1) for i=1:m];
  U2 = [Variable(n_i,1) for i=1:m,j=1:m];

  # Trajectories
  X1 = [Variable(n_s,1) for i=1:m];
  X2 = [Variable(n_s,1) for i=1:m, j = 1:m];
  X3 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m];

  # Epigraph variables
  gamma1 = Variable();
  gamma2 = [Variable() for i=1:m,j=1:m,k=1:m];

  # Problem constraints
  cstr = []
  # System dynamics
  for i1 = 1:m
    cstr += X1[i1] == A_scen[i1]*X_init + B_scen[i1]*U;
    for i2 = 1:m
      cstr += X2[i1,i2] == A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
      for i3 = 1:m
        cstr += X3[i1,i2,i3] == A_scen[i3]*X2[i1,i2] + B_scen[i3]*U2[i1,i2];

        # Final action constraint (a byproduct of using epigraph to minimize)
        cstr += hvcat((2,2),gamma2[i1,i2,i3],X3[i1,i2,i3]',X3[i1,i2,i3],Q[i3]) >= 0
      end
    end
  end

  
  # Cost epigraph constraints
  expected_cost = quad_form(U,R_c);
  
  # Move these out at some point.
  r0 = 0.25;
  max_risk = 0.25;
  alpha = 0.25;

  r1 = [Variable() for i=1:m]
  r2 = [Variable() for i=1:m,j=1:m]

  

  MS = 0;  

  for i1 = 1:m
    # Cost of first action
    c_ijk = quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c);
    MS2 = 0;
    for i2 = 1:m
      c_ijk += quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c);
      MS3 = 0;
      for i3 = 1:m
         c_ijk += gamma2[i1,i2,i3];
         expected_cost += c_ijk*p[i1]*p[i2]*p[i3];
         
         MS3 += p[i3]*max(0,-r2[i1,i2] + gamma2[i1,i2,i3]);
         
         # Condition 2 - check Q??
	 cstr += hvcat((2,2),gamma2[i1,i2,i3], X3[i1,i2,i3]',X3[i1,i2,i3],Q[i3]) >= 0;
      end
      MS3 = r2[i1,i2] + MS3/alpha;
      MS2 = p[i2]*max(0,-r1[i1] + quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c) + MS3);
    end
    MS2 = r1[i1] + MS2/alpha;
    MS = p[i1]*max(0,-r0 + quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c) + MS2);
  end
  MS = r0 + MS/alpha;

  cstr += gamma1 == quad_form(U,R_c) + MS; #Note that this does not include the initial cost X_init'*Q_C*X_init
  cstr += expected_cost <= max_risk;

  # Time to solve the problem!
  problem = minimize(gamma1,cstr);
  solve!(problem, SCSSolver(verbose=false))

  return U, expected_cost, problem 
end


function MPC_look_4(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init, Q)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs

  # Control inputs
  U  = Variable(n_i,1);
  U1 = [Variable(n_i,1) for i=1:m];
  U2 = [Variable(n_i,1) for i=1:m,j=1:m];
  U3 = [Variable(n_i,1) for i=1:m,j=1:m,k=1:m];

  # Trajectories
  X1 = [Variable(n_s,1) for i=1:m];
  X2 = [Variable(n_s,1) for i=1:m, j = 1:m];
  X3 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m];
  X4 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m,ii=1:m];

  # Epigraph variables
  gamma1 = Variable();
  gamma2 = [Variable() for i=1:m,j=1:m,k=1:m,ii=1:m];

  # Problem constraints
  cstr = []
  # System dynamics
  for i1 = 1:m
    cstr += X1[i1] == A_scen[i1]*X_init + B_scen[i1]*U;
    for i2 = 1:m
      cstr += X2[i1,i2] == A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
      for i3 = 1:m
        cstr += X3[i1,i2,i3] == A_scen[i3]*X2[i1,i2] + B_scen[i3]*U2[i1,i2];
        for i4 = 1:m
           cstr += X4[i1,i2,i3,i4] == A_scen[i4]*X3[i1,i2,i3] + B_scen[i4]*U3[i1,i2,i3];
           # Final action constraint (a byproduct of using epigraph to minimize)
           cstr += hvcat((2,2),gamma2[i1,i2,i3,i4],X4[i1,i2,i3,i4]',X4[i1,i2,i3,i4],Q[i4]) >= 0
        end
      end
    end
  end

  
  # Cost epigraph constraints 
  expected_cost = quad_form(U,R_c);
  
  # Move these out at some point.
  r0 = 0.25;
  max_risk = 0.25;
  alpha = 0.25;

  r1 = [Variable() for i=1:m]
  r2 = [Variable() for i=1:m,j=1:m]
  r3 = [Variable() for i=1:m,j=1:m,k=1:m]

  MS = 0;  

  for i1 = 1:m
    # Cost of first action
    c_ijk = quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c);
    MS2 = 0;
    for i2 = 1:m
      c_ijk += quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c);
      MS3 = 0;
      for i3 = 1:m
         c_ijk += quad_form(X3[i1,i2,i3],Q_c) + quad_form(U3[i1,i2,i3],R_c);
         MS4 = 0;
         for i4 = 1:m
           c_ijk += gamma2[i1,i2,i3,i4];
           expected_cost += c_ijk*p[i1]*p[i2]*p[i3]*p[i4];
         
           MS4 += p[i4]*max(0,-r3[i1,i2,i3] + gamma2[i1,i2,i3,i4]);
         
         # Condition 2 - check Q??
	   cstr += hvcat((2,2),gamma2[i1,i2,i3,i4], X4[i1,i2,i3,i4]',X4[i1,i2,i3,i4],Q[i4]) >= 0;
         end
         MS4 = r3[i1,i2,i3] + MS4/alpha;
         MS3 = p[i3]*max(0,-r2[i1,i2] + quad_form(X3[i1,i2,i3],Q_c) + quad_form(U3[i1,i2,i3],R_c) + MS4);
      end
      MS3 = r2[i1,i2] + MS3/alpha;
      MS2 = p[i2]*max(0,-r1[i1] + quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c) + MS3);
    end
    MS2 = r1[i1] + MS2/alpha;
    MS = p[i1]*max(0,-r0 + quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c) + MS2);
  end
  MS = r0 + MS/alpha;

  cstr += gamma1 == quad_form(U,R_c) + MS; #Note that this does not include the initial cost X_init'*Q_C*X_init
  cstr += expected_cost <= max_risk;

  # Time to solve the problem!
  problem = minimize(gamma1,cstr);
  solve!(problem, SCSSolver(verbose=false))

  return U, problem 
end

function MPC_look_5(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init, Q)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs

  # Control inputs
  U  = Variable(n_i,1);
  U1 = [Variable(n_i,1) for i=1:m];
  U2 = [Variable(n_i,1) for i=1:m,j=1:m];
  U3 = [Variable(n_i,1) for i=1:m,j=1:m,k=1:m];
  U4 = [Variable(n_i,1) for i=1:m,j=1:m,k=1:m,ii=1:m];

  # Trajectories
  X1 = [Variable(n_s,1) for i=1:m];
  X2 = [Variable(n_s,1) for i=1:m, j = 1:m];
  X3 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m];
  X4 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m,ii=1:m];
  X5 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m,ii=1:m,jj=1:m];
  
  # Epigraph variables
  gamma1 = Variable();
  gamma2 = [Variable() for i=1:m,j=1:m,k=1:m,ii=1:m,jj=1:m];

  # Problem constraints
  cstr = []
  # System dynamics
  for i1 = 1:m
    cstr += X1[i1] == A_scen[i1]*X_init + B_scen[i1]*U;
    for i2 = 1:m
      cstr += X2[i1,i2] == A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
      for i3 = 1:m
        cstr += X3[i1,i2,i3] == A_scen[i3]*X2[i1,i2] + B_scen[i3]*U2[i1,i2];
        for i4 = 1:m
           cstr += X4[i1,i2,i3,i4] == A_scen[i4]*X3[i1,i2,i3] + B_scen[i4]*U3[i1,i2,i3];
           for i5 = 1:m
              cstr += X5[i1,i2,i3,i4,i5] == A_scen[i5]*X4[i1,i2,i3,i5] + B_scen[i5]*U4[i1,i2,i3,i4];
              # Final action constraint (a byproduct of using epigraph to minimize)
              cstr += hvcat((2,2),gamma2[i1,i2,i3,i4,i5],X5[i1,i2,i3,i4,i5]',X5[i1,i2,i3,i4,i5],Q[i5]) >= 0
           end
        end
      end
    end
  end

  
  # Cost epigraph constraints 
  expected_cost = quad_form(U,R_c);
  
  # Move these out at some point.
  r0 = 0.25;
  max_risk = 0.25;
  alpha = 0.25;

  r1 = [Variable() for i=1:m]
  r2 = [Variable() for i=1:m,j=1:m]
  r3 = [Variable() for i=1:m,j=1:m,k=1:m]
  r4 = [Variable() for i=1:m,j=1:m,k=1:m,ii=1:m]
  MS = 0;  

  for i1 = 1:m
    # Cost of first action
    c_ijk = quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c);
    MS2 = 0;
    for i2 = 1:m
      c_ijk += quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c);
      MS3 = 0;
      for i3 = 1:m
         c_ijk += quad_form(X3[i1,i2,i3],Q_c) + quad_form(U3[i1,i2,i3],R_c);
         MS4 = 0;
         for i4 = 1:m
           c_ijk += quad_form(X4[i1,i2,i3,i4],Q_c) + quad_form(U4[i1,i2,i3,i4],R_c);
           MS5 = 0;
           for i5 = 1:m
             c_ijk += gamma2[i1,i2,i3,i4,i5];
             expected_cost += c_ijk*p[i1]*p[i2]*p[i3]*p[i4]*p[i5];
         
             MS5 += p[i5]*max(0,-r4[i1,i2,i3,i4] + gamma2[i1,i2,i3,i4,i5]);
         
             # Condition 2 - check Q??
	     cstr += hvcat((2,2),gamma2[i1,i2,i3,i4,i5], X5[i1,i2,i3,i4,i5]',X5[i1,i2,i3,i4,i5],Q[i5]) >= 0;
           end
           MS5 = r4[i1,i2,i3,i4] + MS5/alpha;
           MS4 = p[i4]*max(0,-r3[i1,i2,i3] + quad_form(X4[i1,i2,i3,i4],Q_c) + quad_form(U4[i1,i2,i3,i4],R_c) + MS5);
         end
         MS4 = r3[i1,i2,i3] + MS4/alpha;
         MS3 = p[i3]*max(0,-r2[i1,i2] + quad_form(X3[i1,i2,i3],Q_c) + quad_form(U3[i1,i2,i3],R_c) + MS4);
      end
      MS3 = r2[i1,i2] + MS3/alpha;
      MS2 = p[i2]*max(0,-r1[i1] + quad_form(X2[i1,i2],Q_c) + quad_form(U2[i1,i2],R_c) + MS3);
    end
    MS2 = r1[i1] + MS2/alpha;
    MS = p[i1]*max(0,-r0 + quad_form(X1[i1],Q_c) + quad_form(U1[i1],R_c) + MS2);
  end
  MS = r0 + MS/alpha;

  cstr += gamma1 == quad_form(U,R_c) + MS; #Note that this does not include the initial cost X_init'*Q_C*X_init
  cstr += expected_cost <= max_risk;

  # Time to solve the problem!
  problem = minimize(gamma1,cstr);
  solve!(problem, SCSSolver(verbose=false))

  return U, problem 
end

#function that tests MPC_basic_3()
function MPC_run_3(n,m)

  n_s =n;
  n_i =n;
  xi_pts = 0.25*ones(80,m);
  p = 1/m*ones(1,m);
  R_c = 0.0001*eye(n_i,n_i);
  Q_c = 0.0001*eye(n_s,n_s);

  A_scen = fill(eye(n_s),m)
  B_scen = fill(eye(n_s),m)
  X_init = ones(n_s,1)

  R_inv    = inv(R_c);
  Q_c_sqrt   = sqrtm(Q_c);

  Q,toffline = @timed MPC_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c)

  U, tonline = @timed MPC_look_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)

  println(n,", ",toffline,", ",tonline)
end


#function that tests MPC_basic_3()
function MPC_run_N()

  m = 4;
  n_s =5;
  n_i =5;
  xi_pts = 0.25*ones(80,m);
  p = 1/m*ones(1,m);
  R_c = 0.0001*eye(n_i,n_i);
  Q_c = 0.0001*eye(n_s,n_s);

  A_scen = fill(eye(n_s),m)
  B_scen = fill(eye(n_s),m)
  X_init = ones(n_s,1)

  R_inv    = inv(R_c);
  Q_c_sqrt   = sqrtm(Q_c);

  #Run offline code
  println("Setting up");
  Q,toffline = @timed MPC_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c)
  println("Running tests for N=2:5");
  U, t2 = @timed MPC_look_2(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("2 ",t2)
  U, t3 = @timed MPC_look_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("3 ",t3)
  U, t4 = @timed MPC_look_4(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("4 ",t4)
  U, t5 = @timed MPC_look_5(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("5 ",t5)

  println("Running tests for N=2:5");
  U, t2 = @timed MPC_look_2(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("2 ",t2)
  U, t3 = @timed MPC_look_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("3 ",t3)
  U, t4 = @timed MPC_look_4(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("4 ",t4)
  U, t5 = @timed MPC_look_5(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init,Q)
  println("5 ",t5)

end


#function that tests runtime
function MPC_runtime()
    Q,toffline = @timed MPC_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c)
end

