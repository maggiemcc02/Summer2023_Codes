

# In here I will import everything, set parameters, and run codes
# This will be my main code if you will :)




# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import os


# import files


from automatic_solver import *
from classic_arclength import *
from coupling import *
from deflated_continuation_code import *
from deflation import *
from equidistribute import *
from solution_discovery_code import *
from initial_guess_pool import *
from L2_optimal import *
from physical_solver import *
from system_and_jacobian import *



# set parameters


u_0 = 3/2
u_n = 1/2

N = 100

power = 2.5
alpha = 1

physical_tol = 1e-8
boor_tol = 1e-8

damping = True

#eps_0 = 0.05
#eps = eps_0
#eps_f = 1e-6
#delta_eps = 1e-6

eps_list_1 = np.logspace(-1, -6, 41) # set the range of epsilons
eps_list = [eps_list_1[i] for i in range(1, len(eps_list_1))] # take out the initial
E = 0.1 # set the initial

# set grids


grid = np.linspace(0, 1, N)
mesh = np.array([grid]).T
uni_grid = np.array([grid]).T


# set guesses


guess_list = []

guesses = [initial_guess_1, initial_guess_2, initial_guess_3, initial_guess_4, initial_guess_5]




for j in range(5):

  initial_guess = guesses[j]


  guess = []
  for i in range(len(mesh)):
    if i != 0:
      if i != (len(mesh)-1):
        val = [ initial_guess(mesh[i][0]) ]
        guess.append(val)
  guess = np.array( guess )


  guess_list.append(guess)



# call solution discovery



chosen_sols, chosen_mesh = solution_discovery(mesh, uni_grid, guess_list, M_calc_optimal_L2, \
                                                           E, N, damping, alpha, power, [], [], 0, boor_tol, physical_tol, u_0, u_n)


# save the chosen solutions


# create the file path

my_path = '/home/margaretam/github_output/Output_From_MacLachlans_Computer/Continued_Deflation_Codes/output_plots_test2/test2_try2_'




print('My path is', my_path)

for i in range(len(chosen_sols)):

  str_eps_i = str(E)

  str_eps = str_eps_i.replace('0.', '')

  my_plot = 'Solution_' + str(i+1) + '_at_eps_' + str_eps + '.pdf'

  plt.figure()
  plt.plot(chosen_mesh[i], chosen_sols[i], 'blue')
  plt.plot(chosen_mesh[i], [0 for j in range(len(chosen_mesh[i]))], 'blue', marker = "|")
  plt.title('Solution ' + str(i+1) +  " found at eps = " + str(E))
  plt.xlabel('grid')
  plt.ylabel('u(x) approximation')

 # print('The path to save the file is', os.path.join(my_path, my_plot))
  print("The path to save the file is", my_path + my_plot)
  print()

  # plt.savefig(os.path.join(my_path, my_plot))
  plt.savefig(my_path + my_plot)

  plt.close()

 




# now lets call continued deflation



cont_defl_results = continued_deflation(eps_list,  chosen_sols, chosen_mesh, guess, mesh, uni_grid, M_calc_optimal_L2, boor_tol, physical_tol, N, alpha, power, u_0, u_n, guess_list, damping, my_path)






