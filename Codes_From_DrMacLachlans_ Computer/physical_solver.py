
# Here I will include my physical solver codes

# In particular I have
# line search code
# perturbation codes
# physical solver


# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d





from deflation import *
from system_and_jacobian import *










# The Line Search









def inv_quad_interp_line_search(U_old, delta, a_old, a, System, smile, J_xc, F_xc, \
                            mesh, mesh_space, k, k_prime, u_0, u_n, E, solutions, sol_mesh, power, alpha):
  


  max_step = 1 # max newton step we will take 
  step_tol = 1e-12 # smallest damping parameter we can accept
  problem = False
  perturb = False

 
  for j in range(10): # repeat a max of 10 times
 
 
    # compute a_mid and da

    a_mid = (1/2)*(a + a_old)

    da = a - a_old


    

    if da == 0.0:

      break


    # compute ||F|| at each of these values 


    U_a_old = U_old + a_old * delta

    U_a = U_old + a * delta

    U_a_mid = U_old + a_mid * delta


    # compute smiles

    smiles_old = smiley_face(U_a_old, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power)

    smile_old = smiles_old[0]


    smiles_a = smiley_face(U_a, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power)

    smile_a = smiles_a[0]


    smiles_mid = smiley_face(U_a_mid, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power)

    smile_mid = smiles_mid[0]

    # Compute Associated F values


    F_a_old = smile_old * System( U_a_old , mesh, mesh_space, k, u_0, u_n, E)

    F_a = smile_a * System( U_a , mesh, mesh_space, k, u_0, u_n, E)

    F_a_mid = smile_mid * System( U_a_mid , mesh, mesh_space, k, u_0, u_n, E)

    # Compute the Objective function values

    f_old = (1/2) * (norm(F_a_old))**2

    f_a = (1/2) * (norm(F_a))**2

    f_mid = (1/2) * (norm(F_a_mid))**2


    # Compute the derivatives 


    dF = (3*f_a - 4*f_mid + f_old) / da

    dF_old = (-3 * f_old + 4*f_mid - f_a) / da

    d2F = (dF - dF_old) / da


    # check the dF value 

    if dF == 0.0:

      # print('dF = ', dF, 'so we quit')
      # print()
      break

    if d2F == 0.0:

      # print('Mayday, d2F = ', d2F, 'so we gotta quit')
      # print()
      break

    # do the checks and the updates 

    if d2F > 0 :

      a_new = a - (dF / d2F)

    if d2F < 0 :

      a_new = a + (dF / d2F)


    # do breaking checks 

    if np.isnan(a_new): 

      a = None

      break
    
    if a_new > max_step :

      # a is larger then the step we can take
      # so we take a full newton step

      a = 1

      break

    if a_new < step_tol :

      # a is smaller then the step we can take
      # SO we take the midpoint and continue")

      a_new = (1/2) * (a + a_old)

      # problem = True


    # Now do the update for the next iteration

    a_old = a 

    a = a_new


  # do a check - from Schnabel

  # compute macheps 

  # machine_eps = macheps()

  # # calculate the QR factorization of J(xc)

  # Qc, Rc = np.linalg.qr(J_xc)

  # # calculate the condition number of Rc - they use l1 norm 

  # cond_Rc = np.linalg.cond(Rc, 1)

  # # check if we need to perturb

  # # print('condition number of Rc is ', cond_Rc)
  # # print('1/sqrt(macheps) is ', (machine_eps)**(-1/2) )
  # # print()

  # if cond_Rc > ((machine_eps)**(-1/2)):

  #   problem = True

  #   print('We will need to perturb because')
  #   print("condition number of Rc = ", cond_Rc)
  #   print('and macheps = ', macheps)
  #   print() 


  # a differant check

  J_xc_T = np.transpose(J_xc)

  prod_b = np.matmul(J_xc_T, J_xc)  

  prod_c = np.matmul(J_xc_T, F_xc)

  prod_a = np.transpose(prod_c)

  check_2 = (-1) * np.matmul(np.matmul(prod_a, prod_b), prod_c)

  # if check_2 < 0:
  #   # print("We have a descent direction!")
  #   # print()

  if check_2 > 0:

    print('We dont have a descent direction!')
    print()

    problem = True


  if problem: # print out the direction to see whats otg and compute a perturbed step (Schnabel Ch. 6)


    alpha_list = np.linspace(-1/2, 1, 1000)

    norm_F = []

    for ai in alpha_list:

      U_ai = U_old + ai * delta

      smiles_i = smiley_face(U_ai, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power )

      smile_i = smiles_i[0]

      F_Uai = smile_i * System( U_ai , mesh, mesh_space, k, u_0, u_n, E)

      norm_F_Uai = (1/2) * norm(F_Uai)**2 

      norm_F.append(norm_F_Uai)


        
    # # plot the result (positive)

    # plt.plot(alpha_list[333::], np.log10(norm_F[333::]), 'blue')
    # plt.xlabel('alpha')
    # plt.ylabel('(1/2) * || F(U + alpha * delta) || ** 2 ')
    # plt.title("Plot of Norm of F as a Function of Alpha (0 to 1)")
    # plt.show()
    # print()
    # print()


    # # repeat for a larger scale

    # plt.plot(alpha_list, np.log10(norm_F), 'green')
    # plt.xlabel('alpha')
    # plt.ylabel(' (1/2) * || F(U + alpha * delta) || ** 2  ')
    # plt.title("Plot of Norm of F as a Function of Alpha (-1/2 to 1)")
    # plt.show()
    # print()
    # print()

    # # now focus on up to 0.1

    # plt.plot(alpha_list[333:399], np.log10(norm_F[333:399]), 'hotpink')
    # plt.xlabel('alpha')
    # plt.ylabel('(1/2) * || F(U + alpha * delta) || ** 2 ')
    # plt.title("Plot of Norm of F as a Function of Alpha (0, 0.1)")
    # plt.show()
    # print()
    # print()

    # # now focus on up to 0.3

    # plt.plot(alpha_list[333:533], np.log10(norm_F[333:533]), 'purple')
    # plt.xlabel('alpha')
    # plt.ylabel(' (1/2) * || F(U + alpha * delta) || ** 2')
    # plt.title("Plot of Norm of F as a Function of Alpha (0, 0.3)")
    # plt.show()
    # print()
    # print()


    # print("Calling the Pertubation Code")
    # print()

    perturb_step = perturb_system(J_xc, F_xc, smile, System, U_old, mesh, mesh_space,\
                    k, u_0, u_n, E, solutions, sol_mesh, power, alpha)
  

    # print('Damping Parameter for problem plot is :', a)
    # print()

    print("We are going to Perturb the Newton Step")
    print()
    print()
    
   

  else:

    perturb_step = delta 



  return a, perturb_step













# The Perturbation Codes
















# machine epsilon code

def macheps():

  meps = 1 

  loop = True

  iter = 1 

  while loop and iter < 200:

    meps = meps/2

    check_macheps = 1 + meps

    #if ((check_macheps - 1) < 1e-15):
    if check_macheps == 1 :

      loop = False

    iter += 1

  return meps





# Computation of Perturbed System


def perturb_system(J_xc, F_xc, smile, System, U_old, mesh, mesh_space,\
                   k, u_0, u_n, E, solutions, sol_mesh, power, alpha):
  

  # Dennis and Schnabel Ch. 6
  

  # compute macheps 

  machine_eps = macheps()

  # dimensions

  n_perturb = len(J_xc) 

  # compute Hc - perturbed step according to Dennis and Schnabel

  J_xc_T = np.transpose(J_xc)

  term_a = np.matmul(J_xc_T, J_xc) 

  term_b = (n_perturb * machine_eps)**(1/2)

  term_c = norm(term_a, 1) # one norm 

  term_d = (term_b * term_c) * np.identity(n_perturb)

  Hc = term_a + term_d

  # compute the new Newton step

  rhs = (-1) * np.matmul(J_xc_T, F_xc)

  new_delta = np.linalg.solve(Hc, rhs)

  # # print out the newton direction 

  # alpha_list = np.linspace(-1/2, 1, 1000)

  # norm_Fi = []

  # for ai in alpha_list:

  #   U_i = U_old + ai * new_delta

  #   smiles_ui = smiley_face(U_i, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power)

  #   smile_ui = smiles_ui[0]

  #   F_Ui = smile_ui * System( U_i , mesh, mesh_space, k, u_0, u_n, E, solutions, sol_mesh, power, alpha)

  #   norm_F_Ui = (1/2) * norm(F_Ui)**2

  #   norm_Fi.append(norm_F_Ui)

      
  # # plot the result (positive)

  # plt.plot(alpha_list[333::], np.log10(norm_Fi[333::]), 'blue')
  # plt.xlabel('alpha')
  # plt.ylabel('|| F(U + alpha * delta) || ')
  # plt.title("Plot of Norm of F as a Function of Alpha (0 to 1) for Perturbed Code")
  # plt.show()
  # print()
  # print()


  # # repeat for a larger scale

  # plt.plot(alpha_list, np.log10(norm_Fi), 'green')
  # plt.xlabel('alpha')
  # plt.ylabel('1/2 * || F(U + alpha * delta) ||**2 ')
  # plt.title("Plot of Norm of F as a Function of Alpha (-1/2 to 1) for Perturbed Code")
  # plt.show()
  # print()
  # print()

  # # now focus on up to 0.1

  # plt.plot(alpha_list[333:399], np.log10(norm_Fi[333:399]), 'hotpink')
  # plt.xlabel('alpha')
  # plt.ylabel('(1/2) * || F(U + alpha * delta) || ** 2 ')
  # plt.title("Plot of Norm of F as a Function of Alpha (0, 0.1) for Perturbed Code")
  # plt.show()
  # print()
  # print()


  # # now focus on up to 0.3

  # plt.plot(alpha_list[333:533], np.log10(norm_Fi[333:533]), 'purple')
  # plt.xlabel('alpha')
  # plt.ylabel('(1/2) * || F(U + alpha * delta) || ** 2 ')
  # plt.title("Plot of Norm of F as a Function of Alpha (0, 0.3) for Perturbed Code")
  # plt.show()
  # print()
  # print()


  # call the plotting code

  # perturb_plots(J_xc, F_xc, smile, System, U_old, mesh, mesh_space,\
  #                  k, u_0, u_n, E, solutions, sol_mesh, power, alpha)
  
  # return the new delta

  return new_delta



# The Eigenvalue Plots

def perturb_plots(J_xc, F_xc, smile, System, U_old, mesh, mesh_space,\
                   k, u_0, u_n, E, solutions, sol_mesh, power, alpha):
  


  # eigenvalue plots to check if J is positive definite
  
  

  # list of mu's

  mu_list = np.linspace(1e-2, 1e2, 100)
  

  J_xc_T = np.transpose(J_xc)

  term_a = np.matmul(J_xc_T, J_xc) 


  min_e = []

  for j in range(len(mu_list)):

    # compute A = J^TJ + mu * I

    A_pert = term_a + mu_list[j] * np.identity(len(term_a))

    # Compute eigenvalues 

    eig_values, vects = np.linalg.eig(A_pert)

    # take min of eig_values 

    min_eigs = min(eig_values)

    # save min eigenvalue

    min_e.append(min_eigs)

  
  # plot 

  # print('List of min eigenvalues :')
  # print(min_e)
  # print()

  # plt.plot(mu_list, min_e, 'cyan')
  # plt.title("Minimum Eigenvalue vs. Pertubation")
  # plt.xlabel("mu")
  # plt.ylabel('Mimimum Eigenvalue')
  # plt.show()
  # print()
  # print()










# PHYSICAL SOLVER









def new_physical_solver(System, Jacobian, grid, k, k_prime, U0, u_0, u_n, E, tol, max_iter, solutions, sol_mesh, alpha, power, physical_iters,  damping = True):
  """
  Approximates the solution to a specific nonlinear system
  Parameters
  ---------
  System : Vector function
      Is the system of nonlinear equations
  
  Jacobian : Jacobian Matrix of the system

  grid : vector of grid points

  k : Known function in ODE

  k_prime : derivative of k

  U0 :  Vector of initial guesses for each unknown. 
  
  u_0, u_n : boundary conditions

  E : epsilon 
      The pertubation of the problem

  tol : number 
        Tolerance

  max_iter: integer
            Max number of iterations

  damping : real number between 0 and 1
            Damping parameter for Newton's Method
  
  solutions : vector of vectors
            prevously found solutions
  
  sol_mesh : vector of vectors
             corresponding meshes for solutions in solutions

  
  power and alpha : real numbers
                    parameters for deflation operator

  Returns
  ------
  An approximation for the solutions of G(x) = 0 where G is a nonlinear system.
  """


  if damping == False:
    print('NOT going to Damp')
    print()


  # Remove the first and last entries of the grid
  # this is because we already have the solution
  # at the boundary values (u_0 and u_1)
  
  mesh = []
  for i in range(len(grid)):
    if i != 0: # if not first entry
      if i != (len(grid)-1): # if not last entry
        mesh.append([ grid[i][0] ]) # then add it to new vector
  mesh = np.array(mesh)
  

  # Calculate the mesh width (h) for each grid point
  # h_i = x_i - x_i-1

  h_list = []
  for i in range(1, len(grid)):
      h = grid[i][0] - grid[i-1][0]
      h_list.append( [h] )
  mesh_space = np.array(h_list)

  Un = U0 # sets initial guesses as xn for the first iteration of the loop


  
  for n in range(1, max_iter):


    # find smiley face 

    smiles = smiley_face(Un, mesh, mesh_space, u_0, u_n, solutions, sol_mesh, alpha, power)

    smile = smiles[0]
    smile_for_grumpy = smiles[1]

    # find grumpy face 

    grumpy = grumpy_face(Un, mesh, mesh_space, smile_for_grumpy, solutions, sol_mesh, alpha, power)


    # the System

    b = (-1) * smile * System(Un, mesh, mesh_space, k, u_0, u_n, E)

    # The first part of the Jacobian 

    J_term1 = smile * Jacobian(Un, mesh, mesh_space, k, k_prime, u_0, u_n, E)

    # the second part of the jacobian - the outer product 

    J_term2 = np.matmul(System(Un, mesh, mesh_space, k, u_0, u_n, E), grumpy)

    # the Jacobian is term1 + term2

    A = J_term1 + J_term2

    # print out the condition number of A

    # cond_J = np.linalg.cond(A)

    # print('The condition number of J is', cond_J)
    # print()

    # Solve linear system for Newton Step


    delta = np.linalg.solve(A, b)


    # Update : Line Search Damping if asked


    if damping:

      
      damp, new_step = inv_quad_interp_line_search(Un, delta, 0, 1, System, smile, A, -b, \
                      mesh, mesh_space, k, k_prime, u_0, u_n, E, solutions,\
                        sol_mesh, power, alpha)

        
      
      # print('damping parameter is ', damp)
      # print()
      # print('First entry of Newton step : ', delta[0][0])
      # # print('First entry of perturbed step: ', new_step[0][0])
      # print()
      # print('_'*100)
      # print()
      # print()


      U_new = Un + damp * new_step


    # Newton Update

    else:
  
      U_new = Un + delta
    

    if norm(U_new - Un) < tol:
      
      
      # U_new is our solution
      # add in the boundary conditions

      
      solutions_list = [[u_0]]
      for i in U_new:
        solutions_list.append(i)
      solutions_list.append([u_n])

      sol = np.array(solutions_list)



      physical_iters.append(n)




      return sol




    Un = U_new

  print('Max iterations reached and/or sequence is diverging')
  print()
  done = "string"
  return done








