

# In this file I have my solution discovery code
# This code makes use of all initial guesses in my initial guess pool
# to discover solutions

# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



# imports

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
from coupling import *




def solution_discovery(mesh, uni_grid, guess_list, monitor_func, \
                                                           eps, N, damping, alpha, power, found_sols, found_mesh, call_num, boor_tol, physical_tol, u_0, u_n): 
  
  
  results = []


  global solutions
  global sol_mesh

  solutions = []
  sol_mesh = []

  attempt = [found_sols, found_mesh]
                                                             

  print("In solution discovery the initial length of inputted list is", len(found_sols))
  print()                                                           



  for m in range(len(guess_list)):


    solutions = attempt[0]
    sol_mesh = attempt[1]


    # damping = True

    print()
    print('_'*100)
    print('_'*100)
    print('TRYING INITIAL GUESS ' + str(m+1))
    print('_'*100)
    print('_'*100)
    print()

    guess = guess_list[m]


    attempt = ultimate_coupling_unigrid(mesh, uni_grid, guess, monitor_func, eps, N, damping, alpha, power, solutions, sol_mesh, boor_tol, physical_tol, u_0, u_n)

    results.append([len(attempt[0]), attempt[0], attempt[1]]) # number of solutions, solutions, meshes

                                                             



  # return that chosen solution set
                                                             

  print('At the end of solution discovery, the length of our inputted list is', len(found_sols))
  print()
                                                             
  return results[-1][-2], results[-1][-1]



