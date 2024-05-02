# In this file I have the coupling of
# the automatic solver and deflation so that
# we can find multiple solutions
# physical solver has no perturbation codes


# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



from no_perturb_automatic_solver import * # changed
from deflation import *
from no_perturb_physical_solver import * # changed
from equidistribute import *
from system_and_jacobian import *
from classic_arclength import *
from optimal_arc import *
from L2_optimal import *
from numerical_integration import *



def ultimate_coupling_unigrid( grid, uni_grid, guess, mesh_density_calc, E, N,  damping, alpha, power, solutions, sol_mesh, boor_tol, physical_tol, u_0, u_n):


  loop = True

  iter_count = 1


  # start the loop

  while loop: # find 3 solutions

    print()
    print()
    print('~'*100)
    print("DEFLATION NUMBER ", iter_count, '-> We are deflating', len(solutions), "Solutions")
    print('~'*100)
    print()
    print()


    sol, new_mesh = automatic_solver_interpU(F, JF, grid, uni_grid, guess, mesh_density_calc, solutions, sol_mesh, boor_tol, physical_tol, E, N, alpha, power, u_0, u_n)

    if type(sol) == str:
      loop = False
      break


    # check the solution 


    mesh1 = []
    sol1 = []
    for i in range(len(new_mesh)):
      if i != 0: # if not first entry
        if i != (len(new_mesh)-1): # if not last entry
          mesh1.append([new_mesh[i][0] ]) # then add it to new vector
          sol1.append([sol[i][0]])
    mesh1 = np.array(mesh1)
    sol1 = np.array(sol1)
      

    # Calculate the mesh width (h) for each grid point
    # h_i = x_i - x_i-1

    h_list1 = []
    for i in range(1, len(new_mesh)):
        h = new_mesh[i][0] - new_mesh[i-1][0]
        h_list1.append( [h] )
    mesh_space_1 = np.array(h_list1)


    check1 = F(sol1, mesh1, mesh_space_1, k, u_0, u_n, E)


    print('CHECKING THE SOLUTION :', norm(check1))
    print()
    print() 

    solutions.append(sol)
    sol_mesh.append(new_mesh)

    iter_count += 1

    # plt.plot(new_mesh, sol, 'red')
    # plt.plot(new_mesh, [0 for i in range(len(new_mesh))], marker = "|", color = 'darkred')
    # plt.xlabel('Grid Points')
    # plt.ylabel('Solution Function Approximation')
    # plt.title('Final Solution Approx')
    # plt.show()
    # print()
    # print()

  
  return solutions, sol_mesh