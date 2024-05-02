
# Here I will include my physical solver codes

# In particular I have
# line search code
# NOOO perturbation codes
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

  return a, delta





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