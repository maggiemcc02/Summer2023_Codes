

# In this file I have my continued deflation code



# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d

# import files

from automatic_solver import *
from deflation import *
from physical_solver import *
from equidistribute import *
from system_and_jacobian import *
from classic_arclength import *
from optimal_arc import *
from L2_optimal import *
from numerical_integration import *
from initial_guess_pool import *
from solution_discovery_code import *
from coupling import *






def continued_deflation(eps_list, sol_eps, mesh_eps, guess, grid, uni_grid, monitor_func, boor_tol, physical_tol, N, alpha, power, u_0, u_n, guess_list, damping, my_path):




  # list to track number of solutions

  results = []

  call_num = 1

  num_sols = len(sol_eps)

  # set guesses

  guesses = []

  for j in range(num_sols):

    guess = []
    for l in range(len(sol_eps[0])):
      if l != 0:
        if l != (len(sol_eps[0])-1):
          val = [sol_eps[j][l][0]]

          guess.append(val)

    guess = np.array(guess)

    guesses.append(guess)


  yikes_counter = 0
  


  

  for eps in eps_list:

    call_num += 1

    yikes_counter += 1

    if yikes_counter > 100 :

      print("Quitting because yikes occured")
      print()

      return None


    print()
    print()
    print()
    print()
    print('EPSILON = ', eps)
    print()
    print()
    print()
    print()

    global solutions
    global sol_mesh

    solutions = []
    sol_mesh = []




    
    new_sol_eps = []
    new_mesh_eps = []

    

    for i in range(num_sols):

      print()
      print()
      print('-_'*100)
      print('FINDING SOLUTION NUMBER : ', i, 'at eps = ', eps)
      print('_-'*100)
      print()
      print()


      # find a new solution

      new_sol, new_mesh = automatic_solver_interpU(F, JF, mesh_eps[i], uni_grid,\
                                                  guesses[i], monitor_func, solutions,\
                                                  sol_mesh, boor_tol, physical_tol,\
                                                   eps, N, alpha, power, u_0, u_n)
    
      

      if type(new_sol) == str:

        print('We could not find a solution so we move on')
        print()
        continue
      

      # save to list 

      new_sol_eps.append(new_sol)
      new_mesh_eps.append(new_mesh)




    # deflate the solutions out and try and find more

    print()
    print()
    print()
    print('NOW GOING TO DEFLATE')
    print()
    print()
    print()




    new_sol_deflate, new_mesh_deflate = solution_discovery(grid, uni_grid, guess_list, monitor_func, \
                                                           eps, N, damping, alpha, power, new_sol_eps, new_mesh_eps, call_num, boor_tol, physical_tol, u_0, u_n)
    


    num_sols = len(new_sol_eps)

    print()
    print()
    print("WE FOUND", len(new_sol_eps), 'SOLUTIONS')
    print("AFTER DEFLATING WE FOUND", len(new_sol_deflate), "SOLUTIONS")
    print()
    print()

    results.append([eps, len(new_sol_eps), new_sol_eps, new_mesh_eps])


    # create the file path

    #my_path = '/home/margaretam/github_output/Output_From_MacLachlans_Computer/Continued_Deflation_Codes'

    print('Saving plots')
    print()

    for i in range(len(new_sol_eps)):

      str_eps_i = str(eps)

      str_eps = str_eps_i.replace('0.', '')

      my_plot = 'Solution_' + str(i+1) + '_at_eps_' + str_eps + '.pdf'

      plt.figure() 
      plt.plot(new_mesh_eps[i], new_sol_eps[i], 'blue')
      plt.plot(new_mesh_eps[i], [0 for j in range(len(new_mesh_eps[i]))], 'blue', marker = "|")
      plt.title('Solution ' + str(i+1) +  " found at eps = " + str(eps))
      plt.xlabel('grid')
      plt.ylabel('u(x) approximation')

      # plt.savefig(os.path.join(my_path, my_plot))
      plt.savefig(my_path + my_plot)
      plt.close()




    # update the guesses

    guesses = []

    for j in range(num_sols):
     
      guess = []
     
      for l in range(len(new_sol_eps[0])):
        
        if l != (len(new_sol_eps[0])-1):
          
          if l != 0:

            val = [new_sol_eps[j][l][0]]
            
            guess.append(val)



      guess = np.array(guess)

      guesses.append(guess)

      

    # update the guesses

    mesh_eps = []
    for j in range(len(new_mesh_eps)):
      mesh_eps.append(new_mesh_eps[j])



    # update epsilon

    #eps -= delta_eps


  # print the numbers of solutions

  print()
  print()
  print('_'*100)
  print()

  for i in range(len(results)):

    print('For eps = ', results[i][0] , 'we found', results[i][1], 'solutions')
    print()

  print()
  print('_'*100)
  print()
  print()




  

  return results











































