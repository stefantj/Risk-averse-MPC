using Convex
using SCS
using Debug



#Function that computes control for look-ahead value of 3
function MPC_basic_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init)

  # Infer problem size
  m   =  size(B_scen,1);    # Number of realizations
  n_s =  size(B_scen[1],1); # Number of state variables
  n_i =  size(B_scen[1],2); # Number of control inputs


  # Time to define Optimization variables
  
  # These are Decision variables for the LMI that describes optimality
  Q = [Semidefinite(n_s) for i=1:m];
  G = [Variable(n_s,n_s) for i=1:m];
  Y = [Variable(n_i,n_s) for i=1:m];

  # These are Control inputs
  U  = Variable(n_i,1);
  U1 = [Variable(n_i,1) for i=1:m];
  U2 = [Variable(n_i,1) for i=1:m,j=1:m];

  # These are trajectories
  X1 = [Variable(n_s,1) for i=1:m];
  X2 = [Variable(n_s,1) for i=1:m, j = 1:m];
  X3 = [Variable(n_s,1) for i=1:m, j = 1:m, k = 1:m];

  # These are epigraph variables
  gamma1 = Variable();
  gamma2 = [Variable() for i=1:m,j=1:m,k=1:m];


  # Define some useful constants that are used a lot
  R_inv    = inv(R_c);
  R_c_sqrt   = sqrtm(R_c);
  Q_c_sqrt   = sqrtm(Q_c);


  # Define the problem constraints

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

  # System dynamics
  for i1 = 1:m
    cstr += X1[i1] == A_scen[i1]*X_init + B_scen[i1]*U;
    for i2 = 1:m
      cstr += X2[i1,i2] == A_scen[i2]*X1[i1] + B_scen[i2]*U1[i1];
      for i3 = 1:m
        cstr += X3[i1,i2,i3] == A_scen[i3]*X2[i1,i2] + B_scen[i2]*U2[i1,i2];

        # Final action constraint (a byproduct of using epigraph to minimize)
        cstr += hvcat((2,2),gamma2[i1,i2,i3],X3[i1,i2,i3]',X3[i1,i2,i3],Q[i3]) >= 0
      end
    end
  end

  
  # Cost epigraph constraints 
  for i1 = 1:m
    for i2 = 1:m
      z = 0;
      for i3 = 1:m
         z += quad_form(X3[i1,i2,i3],p[i3]*Q_c) + quad_form(U2[i1,i2],p[i3]*R_c) +gamma2[i1,i2,i3]
	cstr += hvcat((2,2),gamma2[i1,i2,i3], X3[i1,i2,i3]',X3[i1,i2,i3],Q_c) >= 0;
      end
      # Enforce epigraph
      cstr += gamma1 >= z
    end
  end


  # Time to solve the problem!
  println("Parsing problem");
  @time problem = minimize(gamma1,cstr);
  println("Solving problem");
  @time solve!(problem, SCSSolver(verbose=false))


  println("Done with basic MPC step");
  return U,Q,G,Y,gamma2,X3 
end




#function that tests MPC_basic_3()
@debug function MPC_test_3()

  
  m = 4;
  n_s =3;
  n_i =3;
  xi_pts = 0.25*ones(80,m);
  p = 1/m*ones(1,m);
  R_c = 0.0001*eye(n_i,n_i);
  Q_c = 0.0001*eye(n_s,n_s);

  A_scen = fill(eye(n_s),m)
  B_scen = fill(eye(n_s),m)
  X_init = ones(n_s,1)

  R_inv    = inv(R_c);
  Q_c_sqrt   = sqrtm(Q_c);

  U,Q,G,Y,gamma2,X3 = @time MPC_basic_3(A_scen,B_scen,xi_pts,p,R_c,Q_c,X_init)

  println("Checking that Q are PD")
  for i = 1:m
    tmp = eigmin(Q[i].value);
    if(tmp <= -1e-6)
      print_with_color(:red, string("Q[",i,"] not PD: eigmin = ",tmp,"\n"));
    end
  end

  println("Checking that LMI (17) is PSD");
  for i = 1:m
      Q_diag = kron(Q[i].value,eye(m))

      b12 = zeros(m*n_s,n_i);
      b13 = zeros(m*n_s,n_s);
      b22 = R_inv;
      b23 = zeros(n_i,n_s);
      b24 = -(Y[i].value)
      b33 = eye(n_s);
      
      #sigma matrix not necessary, since for our example everything is constant
      b14 = -0.5*(A_scen[1]*G[1].value + B_scen[1]*Y[1].value);
      for j = 2:m
        b14 = vcat(b14, -0.5*(A_scen[j]*G[j].value + B_scen[j]*Y[j].value));
      end

      b34 = -Q_c_sqrt*G[i].value;
      b44 = -Q[i].value + G[i].value + G[i].value';

      k= hvcat((4,4,4,4), Q_diag, b12, b13, b14,
                           b12', b22, b23, b24,
			   b13', b23', b33, b34,
			   b14', b24', b34', b44)

      tmp = eigmin(k);
      if(tmp < -1e-6)
        print_with_color(:red,string("LMI ",i," is not PSD: eigmin = ",tmp,"\n"))
      end
  end

  println("Checking that condition 3 holds")
  for i1 = 1:m
    for i2 = 1:m
      for i3 = 1:m
        tmp = hvcat((2,2),gamma2[i1,i2,i3].value,(X3[i1,i2,i3].value)',(X3[i1,i2,i3].value),(Q[i3].value));
	tmp_e = eigmin(tmp);
	if(tmp_e < -1e-6)
	  print_with_color(:red, string("C3 violated for (",i1,",",i2,",",i3,"): eigmin = ", tmp_e,"\n"));
	end
      end
    end
  end


  println("Checks complete!")

end
