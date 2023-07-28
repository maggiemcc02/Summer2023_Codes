

# In this file I have my automatic solver
# The automatic solver is the coupling of the physical solver and de Boor

# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



from automatic_solver import *
from deflation import *
from physical_solver import *
from equidistribute import *
from system_and_jacobian import *
from classic_arclength import *
from optimal_arc import *
from L2_optimal import *
from numerical_integration import *



def automatic_solver_interpU(System, Jacobian, mesh, uni_grid, guess, mesh_density_calc, solutions, sol_mesh, boor_tol, physical_tol, E, N, alpha, power ,u_0, u_n, finding_sol3 = False, damping = True, damp_guess = 0.5):


  """
  System : Nonlinear system we are solving

  Jacobian : Jacobian of system 

  mesh : mesh we are solving problem on

  uni_grid : uniform mesh for interval we are solving problem on (xi)

  guess : veector of initial guesses for physical solver

  solutions : known solutions (if exists)

  sol_mesh : mesh r was found on

  tol : tolerance for loop - stop when differance between meshes are less than tol

  E : epsilon for the problem

  n : number of grid points

  damping : damping parameter

  power and alpha : real numbers 
                    parameters for the deflation operator

  """

  physical_iters = []

  # print()
  # print('_'*50)
  # print('Initial Physical Solve')
  # print('_'*50)
  # print()

  U = new_physical_solver(System, Jacobian, mesh, k, k_prime, guess, u_0, u_n, E, physical_tol, 100, solutions, sol_mesh, alpha, power, physical_iters)

  # set new initial guess as U

  if type(U) == str:

    return "quit", []

  guess = []
  for i in range(len(U)):
    if i != 0:
      if i != (len(U)-1):
        val = [  U[i][0]  ]
        guess.append(val)
  guess = np.array( guess )


  loop = True

  # loop to couple codes

  global iter
  iter = 1


  while (loop and (iter < 100)): 


    # update iter counter
    
    iter += 1

    # compute mesh density function

    M_vals = mesh_density_calc(U, mesh)

    # equidistribute (mesh solver with de Boor)


    new_mesh = equidistribute(mesh, M_vals, uni_grid, 0, 1, N)

    # compute change in mesh

    check = norm(new_mesh - mesh) 
    
    # This is a dangerous check because 
    # the mesh could change slightly while the solution changes greatly

    if check < boor_tol:
      loop = False
      break


    # interpolate U onto the mesh 

    # make U and its mesh 1D

    U_1D = []
    m_1D = []

    for j in range(len(U)):

      U_1D.append(U[j][0])
      m_1D.append(mesh[j][0])

    U_interp = np.interp(new_mesh, m_1D, U_1D)


    # update the mesh

    mesh = new_mesh

    # set current U as guess
  
    guess = []
    for i in range(len(U_interp)):
      if i != 0:
        if i != (len(U_interp)-1):
          val = [  U_interp[i][0]  ]
          guess.append(val)
    guess = np.array( guess )


    new_U = new_physical_solver( System, Jacobian, mesh, k, k_prime, guess,\
                                 u_0, u_n, E, physical_tol, 100, solutions, sol_mesh, alpha, power, physical_iters)

    # set new_U as U

    U = new_U

    if type(new_U) == str:

      print("The physical solver failed so we are 'out' of solutions")
      return "yikes", []


  # check if we diverged

  if iter >= 100:

    print('We have diverged in the physical + de Boor loop')
    print()
    print()
    return 'yikes', []

  # create final grid 

  #final_mesh = []
  #for i in range(len(mesh)):
    #final_mesh.append(mesh[i][0])
  #final_mesh = np.array([final_mesh]).T

  final_mesh = mesh

  # how many iterations?

  print()
  print()
  print('This solution required ', iter, 'de Boor + Physical Solve Loops')
  print()
  print('The physical solver iteration counts are : ')
  print(physical_iters)
  print()
  print()


  return U, final_mesh






