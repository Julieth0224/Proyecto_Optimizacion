prob = optimproblem; 
prob.Objective = 200*x11 + 300*x21 + 500*x31 + 200*x12 + 300*x22 + 500*x32 
                + 200*x13 + 300*x23 + 500*x33 + 200*x14 + 300*x24 + 500*x34 
                + 200*x15 + 300*x25 + 500*x35 + 200*x16 + 300*x26 + 500*x36 
                + 200*x17 + 300*x27 + 500*x37 
                -((3*X1+1)(X1-5) + (3*X2+1)(X2-5) + (3*X3+1)(X3-5) + (3*X4+1)(X4-5) + (3*X5+1)(X5-5) + (3*X6+1)(X6-5) + (3*X7+1)(X7-5))
                -(20*(X1+5) + 20*(X2+5) + 20*(X3+5) + 20*(X4+5) + 20*(X5+5) + 20*(X6+5) + 20*(X7+5)); 
prob.Constraints.cons1 = xij <= 6;
prob.Constraints.cons2 = Xj - 5 >= Lj;