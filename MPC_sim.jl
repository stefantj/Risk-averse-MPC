# MPC_sim.jl
#
# This file contains functions for testing risk-averse MPC
# Created by Stefan Jorgensen for AA203
# June 6, 2015
#

include("MPC_control.jl");
using Debug

# Store data in csv:
# [ Status                     ] [ State/action                               ] [ Performance         ]
# iteration, feasible, rand_ind, X(1), ... , X(n_state), u(1), ... ,u(n_input), r, c_i, gamma_i, E(c)
function write_data(filename, iter, feas, realization, X, U, risk, cost, gamma, expect_cost)
  file_handle = open(filename, "a");
  write(file_handle, string(iter, ", ", feas, ", ", realization, ", "));
  for i = 1:length(X)
    write(file_handle, string(X[i],", "));
  end
  for i = 1:length(U)
    write(file_handle, string(U[i], ", "));
  end
  write(file_handle, string(risk, ", ", cost, ", ", gamma, ", ", expect_cost, "\n"));
  close(file_handle);
end



# MPC_sim() is a simple way to run the example from 
# the paper by Chow, Pavone
# 
@debug function MPC_sim()
  # Start by assuming that we're in the 2d case and use the example from the paper
  # Simulation length
  num_trials = 2;
  max_iters = 50; 

  # Define the system:
  m = 4;
  n_input = 2;
  n_state = 2;
  
  # Expected risk
  p = 1/m*ones(1,m);
  # Quadratic cost functions
  R_c = 0.0001*eye(n_input)
  Q_c = 0.25*eye(n_state);

  # The example way of generating random outcomes doesn't really make sense to me
  # I'm ignoring it and doing what is intuitive.
  p_cdf = cumsum(p');
  r_tmp = rand(max_iters*num_trials,1);
  rand_indices = zeros(max_iters*num_trials,1);
  for iter = 1:m
    rand_indices += (r_tmp .< p_cdf[iter])
  end

  #TODO: This segment is a mess.

  # Risk levels:
#  alpha = linspace(0,1,num_trials);
  alpha = 0.25;

  # Enumerate vertices in the polytope
  # (I have not verified the correctness of these)
  A_Cvar = [ [ones(1,m) zeros(1,m)];
            -[ones(1,m) zeros(1,m)];
	     [eye(m) diagm(p[:])*(ones(m,1)*p-eye(m))];
	    -[eye(m) diagm(p[:])*(ones(m,1)*p-eye(m))];
             [zeros(m,m) eye(m)];
	    -[zeros(m,m) eye(m)]; ];

  b_Cvar = [1;
            -0.9999; 
	         p'; 
	    -0.999*p';
	 alpha*ones(m,1);
	    zeros(m,1);];

  #TODO: I do not have a reasonable replacement for this code:
  # P = polytope(A_Cvar, b_Cvar)
  # ext = extreme(P);
  # ext = ext(:,1:m)
  # ind = find(abs(sum(ext,2)-1)<1e-3);
  # xi_pts = ext(ind,:)
  # xi=xi_pts(1,:);
  xi_pts = 0.25*ones(40,m);

  # Dynamics:
  A_scen = fill(eye(n_state), m);
  B_scen = fill(zeros(n_state,n_input),m);

  A_scen[1] = [2    0.5;
               -0.5   2];

  A_scen[2] = [-0.1564 -0.0504;
               -0.0504 -0.1904];

  A_scen[3] = [1.50000 -0.3000;
               0.20000 1.50000];

  A_scen[4] = [0.57680 0.28590;
               0.28590 0.74990];
  B_scen[1] = [3 0;
               0 2];

  B_scen[2] = [-0.9540 0;
               -0.7773 0.1852];
  B_scen[3] = [2 0;
               0 2];
  B_scen[4] = [-0.2587 -0.9364;
               0.4721  0];

  # Initial values
  X_init = ones(n_state,1);
  r_k = 50;
  curr_cost_obj = X_init'*Q_c*X_init;

  # State trackers
  X_traj = fill(X_init, num_trials, max_iters);
  r_traj = fill(0.0, num_trials, max_iters);
  X_total = fill(X_init, num_trials);
  r_total = fill(r_k, num_trials);


  # Store data in csv:
  # [ Status                     ] [ State/action                               ] [ Performance         ]
  # iteration, feasible, rand_ind, X(1), ... , X(n_state), u(1), ... ,u(n_input), r, c_i, gamma_i, E(c)
  f_data = open("data.csv","w");
  close(f_data);

  # Solve offline version of the problem
  # if infeasible, then don't continue

  for trials = 1:num_trials
    # Solve off-line version of the problem:
    Q_offline, pr = MPC_offline(A_scen, B_scen, xi_pts, p, R_c, Q_c);
    if(pr.status == :Infeasible)
      println("Infeasible start - unable to solve offline problem.");
      continue
    end

    for iter = 1:(max_iters-1)
      println("Trial ", trials, ", iterate ", iter);
      # For now, only bother with N = 2

      @time ui, Ec,pr = MPC_look_2(A_scen, B_scen, xi_pts, p, R_c, Q_c, X_traj[trials,iter], Q_offline)
      
      feas_status = 0;
      if(pr.status == :Infeasible)
          # Print infeasible
	  println("Infeasible iterate - unable to solve at ", iter);
          feas_status = -1;
      else
	  if(pr.status == :Optimal)
	      # Note that things are good
	      feas_status = 1;
	  else
	      # Note that things may be bad
	      feas_status = 0;
	  end
          # Continue with simulation
          r_i = rand_indices[(trials-1)*max_iters + iter];
          X_traj[trials,iter+1] = A_scen[r_i]*X_traj[trials,iter] + B_scen[r_i]*ui;
          cost_i = X_traj[trials,iter+1]'*Q_c*X_traj[trials,iter+1] + ui'*R_c*ui;

	  # Print data to file
	  write_data("data.csv", iter, feas_status, r_i, X_traj[trials,iter], ui,0, cost_i[1], pr.optval, 0); 
	 
	  if(norm(X_traj[trials,iter+1]) < 1e-4)
             println("Found the origin! (norm = ", norm(X_traj[trials,iter+1]), ")");
	     continue
          end
      end
    end
  end
end





